#include "clr.h"
#include "mixtureClasses.h"
#include "distanceClasses.h"

#include <Eigen/Core>
#include <LBFGS.h>

Rcpp::NumericVector GetMixtureDensity(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::DataFrame &inputModel,
    const bool logScale)
{
  unsigned int numPoints = inputValues.size();
  Rcpp::NumericVector meanValues = inputModel["mean"];
  Rcpp::NumericVector precisionValues = inputModel["precision"];
  Rcpp::NumericVector mixingValues = inputModel["mixing"];
  std::vector<double> workVector;

  Rcpp::NumericVector outputValues(numPoints);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double workScalar = GetLogDensityValue(inputValues[i], meanValues, precisionValues, mixingValues, workVector);

    if (!logScale)
      workScalar = std::exp(workScalar);

    outputValues[i] = workScalar;
  }

  return outputValues;
}

Rcpp::NumericVector GetMean(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::List &inputData,
    const double alpha)
{
  unsigned int numInputs = inputData.size();
  unsigned int numPoints = inputValues.size();
  Rcpp::NumericVector outputValues(numPoints);
  std::vector<double> workVector;
  Rcpp::DataFrame workModel;
  Rcpp::NumericVector meanValues;
  Rcpp::NumericVector precisionValues;
  Rcpp::NumericVector mixingValues;
  std::vector<double> sampleValues(numInputs);
  unsigned int numUselessInputs = std::floor(alpha * numInputs);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double inputValue = inputValues[i];

    for (unsigned int j = 0;j < numInputs;++j)
    {
      workModel = inputData[j];
      meanValues = workModel["mean"];
      precisionValues = workModel["precision"];
      mixingValues = workModel["mixing"];

      sampleValues[j] = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, workVector);
    }

    std::nth_element(sampleValues.begin(), sampleValues.begin() + numUselessInputs, sampleValues.end());

    double outputValue = 0.0;

    for (unsigned int j = numUselessInputs;j < numInputs;++j)
      outputValue += sampleValues[j];

    outputValues[i] = std::exp(outputValue / (numInputs - numUselessInputs));
  }

  return outputValues;
}

Rcpp::NumericVector GetSquaredDistancesToMean(
    const Rcpp::List &inputData,
    const Rcpp::LogicalVector &subsetValues,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numInputs = inputData.size();
  unsigned int numPoints = nodeValues.size();
  Rcpp::NumericVector squaredDistancesToMean(numInputs);

  unsigned int numSubset = 0;
  for (unsigned int i = 0;i < numInputs;++i)
  {
    if (subsetValues[i])
      ++numSubset;
  }

  std::vector<double> referenceMeanValues, referencePrecisionValues, referenceMixingValues;
  unsigned int numComponents = BuildReferenceModel(inputData, referenceMeanValues, referencePrecisionValues, referenceMixingValues);

  Rcpp::DataFrame workModel;
  Rcpp::NumericVector workMeanValues, workPrecisionValues, workMixingValues;
  std::vector<double> workVector;

  // Get distances to mean
  for (unsigned int i = 0;i < numInputs;++i)
  {
    double totalLogIntegral = 0.0;
    double totalSquareLogIntegral = 0.0;

    for (unsigned int j = 0;j < numComponents;++j)
    {
      double logIntegral = 0.0;
      double squareLogIntegral = 0.0;
      double referenceMixingValue = referenceMixingValues[j];

      for (unsigned int k = 0;k < numPoints;++k)
      {
        double inputValue = referenceMeanValues[j] + std::sqrt(2.0 / referencePrecisionValues[j]) * nodeValues[k];
        double weightValue = weightValues[k];

        double logValue = 0.0;

        for (unsigned int l = 0;l < numInputs;++l)
        {
          if (!subsetValues[l])
            continue;

          workModel = inputData[l];
          workMeanValues = workModel["mean"];
          workPrecisionValues = workModel["precision"];
          workMixingValues = workModel["mixing"];

          double workValue = GetLogDensityValue(inputValue, workMeanValues, workPrecisionValues, workMixingValues, workVector);

          if (l == i)
            logValue += (numSubset - 1.0) * workValue;
          else
            logValue -= workValue;
        }

        logValue /= numSubset;

        if (!subsetValues[i])
        {
          workModel = inputData[i];
          workMeanValues = workModel["mean"];
          workPrecisionValues = workModel["precision"];
          workMixingValues = workModel["mixing"];
          logValue += GetLogDensityValue(inputValue, workMeanValues, workPrecisionValues, workMixingValues, workVector);
        }

        logIntegral += weightValue * logValue;
        squareLogIntegral += weightValue * logValue * logValue;
      }

      totalLogIntegral += referenceMixingValue * logIntegral;
      totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
    }

    squaredDistancesToMean[i] = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
  }

  return squaredDistancesToMean;
}

Rcpp::NumericVector GetSquaredDistanceMatrix(
    const Rcpp::List &inputData,
    const std::vector<double> &nodeValues,
    const std::vector<double> &weightValues)
{
  Rcpp::DataFrame workTibble;
  GaussianMixture workMixture;
  double firstMeanValue, secondMeanValue;
  unsigned int numInputs = inputData.size();
  Rcpp::NumericVector outputValues(numInputs * (numInputs - 1) / 2);
  SquaredBayesDistance sqDist;
  sqDist.SetQuadraturePoints(nodeValues);
  sqDist.SetQuadratureWeights(weightValues);

  // Set up parameters
  LBFGSpp::LBFGSParam<double> param;
  param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

  // Create solver and function object
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBacktracking> solver(param);

  // Initial guess
  Eigen::VectorXd x = Eigen::VectorXd::Zero(1);
  // x will be overwritten to be the best point found
  double fx;

  for (unsigned int i = 0;i < numInputs - 1;++i)
  {
    if (i != 0)
      continue;

    workTibble = inputData[i];
    workMixture.SetInput(workTibble);
    firstMeanValue = workMixture.GetMean();
    sqDist.SetInput1(workTibble);

    Rcpp::checkUserInterrupt();

    for (unsigned int j = i + 1;j < numInputs;++j)
    {
      // Rcpp::Rcout << i << " " << j << std::endl;
      if (j != 107)
        continue;

      workTibble = inputData[j];
      workMixture.SetInput(workTibble);
      secondMeanValue = workMixture.GetMean();
      sqDist.SetInput2(workTibble);

      x[0] = firstMeanValue - secondMeanValue;

      Rcpp::Rcout << x << std::endl;
      int niter = solver.minimize(sqDist, x, fx);
      Rcpp::Rcout << x << std::endl;

      // bool success = optim::lbfgs(x, ComputeSquaredDistance, &fData, algoSettings);

      // Rcpp::Rcout << x << std::endl;

      // if (!success)
      //   Rcpp::stop("Optimization failed when computing distances.");

      // double fx = ComputeSquaredDistance(x, &grad, &fData);

      outputValues[numInputs * i - (i + 1) * i / 2 + j - i - 1] = fx;
    }
  }

  return outputValues;
}

double GetMeanRawMoment(
    const Rcpp::List &inputData,
    const unsigned int order,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numInputs = inputData.size();
  unsigned int numPoints = nodeValues.size();

  Rcpp::DataFrame inputModel;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector <double> workVector;

  // First define the reference Gaussian mixture wrt which integration will take place with GH.
  // The best for accuracy is to concatenate all input mixture models into one.
  // Note however that it is computationally expensive.

  std::vector<double> referenceMeanValues, referencePrecisionValues, referenceMixingValues;
  unsigned int numComponents = BuildReferenceModel(inputData, referenceMeanValues, referencePrecisionValues, referenceMixingValues);

  // Now integrate to get normalization constant
  double totalIntegralValue = 0.0;
  double kthTotalIntegralValue = 0.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    double referenceMixingValue = referenceMixingValues[i];
    double integralValue = 0.0;
    double kthIntegralValue = 0.0;

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double nodeValue = nodeValues[j];
      double weightValue = weightValues[j];
      double inputValue = referenceMeanValues[i] + std::sqrt(2.0 / referencePrecisionValues[i]) * nodeValue;

      double logarithmArithmeticMean = 0.0;
      double arithmeticMeanLogarithm = 0.0;

      for (unsigned int k = 0;k < numInputs;++k)
      {
        inputModel = inputData[k];
        meanValues = inputModel["mean"];
        precisionValues = inputModel["precision"];
        mixingValues = inputModel["mixing"];

        double logValue = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, workVector);
        logarithmArithmeticMean += logValue;
        arithmeticMeanLogarithm += std::exp(logValue);
      }

      logarithmArithmeticMean /= numInputs;
      arithmeticMeanLogarithm /= numInputs;
      arithmeticMeanLogarithm = std::log(arithmeticMeanLogarithm);

      double workScalar = std::exp(logarithmArithmeticMean - arithmeticMeanLogarithm);

      integralValue += weightValue * workScalar;
      if (order > 0)
        kthIntegralValue += weightValue * workScalar * std::pow(inputValue, (double)order);
    }

    totalIntegralValue += referenceMixingValue * integralValue;
    if (order > 0)
      kthTotalIntegralValue += referenceMixingValue * kthIntegralValue;
  }

  return (order == 0) ? totalIntegralValue : kthTotalIntegralValue / totalIntegralValue;
}

double GetLogDensityValue(
    const double inputValue,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    std::vector<double> &workVector)
{
  unsigned int numComponents = meanValues.size();
  workVector.resize(numComponents);

  // Find component with smallest exponent
  // and compute log densities of components
  unsigned int indexOfMinimalExponent = 0;
  double minimalExponent = 0.0;
  for (unsigned int i = 0;i < numComponents;++i)
  {
    double workValue = precisionValues[i] * (inputValue - meanValues[i]) * (inputValue - meanValues[i]);

    if (workValue < minimalExponent || i == 0)
    {
      minimalExponent = workValue;
      indexOfMinimalExponent = i;
    }

    workVector[i] = 0.5 * std::log(precisionValues[i] / (2.0 * M_PI)) - workValue / 2.0;
  }

  double insideLogValue = 1.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    if (i == indexOfMinimalExponent)
      continue;

    insideLogValue += std::exp(workVector[i] - workVector[indexOfMinimalExponent]) * mixingValues[i] / mixingValues[indexOfMinimalExponent];
  }

  // Compute log of mixture density
  return std::log(mixingValues[indexOfMinimalExponent]) + workVector[indexOfMinimalExponent] + std::log(insideLogValue);
}

double GetSquaredDistance(
    const Rcpp::NumericVector &firstMeanValues,
    const Rcpp::NumericVector &firstPrecisionValues,
    const Rcpp::NumericVector &firstMixingValues,
    const Rcpp::NumericVector &secondMeanValues,
    const Rcpp::NumericVector &secondPrecisionValues,
    const Rcpp::NumericVector &secondMixingValues,
    const std::vector<double> &referenceMeanValues,
    const std::vector<double> &referencePrecisionValues,
    const std::vector<double> &referenceMixingValues,
    const unsigned int numComponents,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues,
    std::vector<double> &workVector)
{
  unsigned int numPoints = nodeValues.size();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    double referenceMixingValue = referenceMixingValues[i];

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double inputValue = referenceMeanValues[i] + std::sqrt(2.0 / referencePrecisionValues[i]) * nodeValues[j];
      double weightValue = weightValues[j];

      double logValue1 = GetLogDensityValue(inputValue, firstMeanValues, firstPrecisionValues, firstMixingValues, workVector);
      double logValue2 = GetLogDensityValue(inputValue, secondMeanValues, secondPrecisionValues, secondMixingValues, workVector);

      logIntegral += weightValue * (logValue1 - logValue2);
      squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
  }

  return totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
}

unsigned int BuildReferenceModel(
    const Rcpp::List &inputData,
    std::vector<double> &referenceMeanValues,
    std::vector<double> &referencePrecisionValues,
    std::vector<double> &referenceMixingValues)
{
  unsigned int numInputs = inputData.size();
  Rcpp::DataFrame inputModel;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;

  referenceMeanValues.clear();
  referencePrecisionValues.clear();
  referenceMixingValues.clear();

  unsigned int numComponents = 0;

  for (unsigned int i = 0;i < numInputs;++i)
  {
    inputModel = inputData[i];
    unsigned int internalNumComponents = inputModel.nrows();

    meanValues = inputModel["mean"];
    referenceMeanValues.insert(referenceMeanValues.end(), meanValues.begin(), meanValues.end());

    precisionValues = inputModel["precision"];
    referencePrecisionValues.insert(referencePrecisionValues.end(), precisionValues.begin(), precisionValues.end());

    mixingValues = inputModel["mixing"];
    for (unsigned int j = 0;j < internalNumComponents;++j)
      mixingValues[j] /= numInputs;
    referenceMixingValues.insert(referenceMixingValues.end(), mixingValues.begin(), mixingValues.end());

    numComponents += internalNumComponents;
  }

  return numComponents;
}
