#include "clr.h"

Rcpp::NumericVector GetMean(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::List &mixtureModels,
    const double alpha)
{
  unsigned int numSamples = mixtureModels.size();
  unsigned int numPoints = inputValues.size();
  Rcpp::NumericVector outputValues(numPoints);
  std::vector<double> logDensityValues;
  Rcpp::DataFrame mixtureModel;
  Rcpp::NumericVector meanValues;
  Rcpp::NumericVector precisionValues;
  Rcpp::NumericVector mixingValues;
  std::vector<double> sampleValues(numSamples);
  unsigned int numUselessSamples = std::floor(alpha * numSamples);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double inputValue = inputValues[i];

    for (unsigned int j = 0;j < numSamples;++j)
    {
      mixtureModel = mixtureModels[j];
      meanValues = mixtureModel["mean"];
      precisionValues = mixtureModel["precision"];
      mixingValues = mixtureModel["mixing"];

      sampleValues[j] = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);
    }

    std::nth_element(sampleValues.begin(), sampleValues.begin() + numUselessSamples, sampleValues.end());

    double outputValue = 0.0;

    for (unsigned int j = numUselessSamples;j < numSamples;++j)
      outputValue += sampleValues[j];

    outputValues[i] = outputValue / (numSamples - numUselessSamples);
  }

  return outputValues;
}

Rcpp::NumericVector GetSquaredDistancesToMean(
    const Rcpp::List &mixtureModels,
    const Rcpp::LogicalVector &subsetValues,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numNodes = nodeValues.size();
  unsigned int numSamples = mixtureModels.size();
  Rcpp::NumericVector squaredDistancesToMean(numSamples);

  unsigned int numSubset = 0;

  for (unsigned int i = 0;i < numSamples;++i)
    if (subsetValues[i])
      ++numSubset;

  Rcpp::DataFrame referenceModel;
  Rcpp::NumericVector referenceMeanValues, referencePrecisionValues,referenceMixingValues;

  Rcpp::DataFrame mixtureModel;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector<double> logDensityValues;

  // First, recompute optimal reference in this subset of data
  unsigned int pos = 0;
  unsigned int referenceIndex = 0;
  double referenceSquaredDistance = 0;

  for (unsigned int i = 0;i < numSamples;++i)
  {
    if (!subsetValues[i])
      continue;

    referenceModel = mixtureModels[i];
    referenceMeanValues = referenceModel["mean"];
    referencePrecisionValues = referenceModel["precision"];
    referenceMixingValues = referenceModel["mixing"];
    unsigned int referenceNumComponents = referenceModel.nrows();

    double totalLogIntegral = 0.0;
    double totalSquareLogIntegral = 0.0;

    for (unsigned int j = 0;j < referenceNumComponents;++j)
    {
      double logIntegral = 0.0;
      double squareLogIntegral = 0.0;
      double referenceMixingValue = referenceMixingValues[j];

      for (unsigned int k = 0;k < numNodes;++k)
      {
        double logValue = 0.0;
        double inputValue = referenceMeanValues[j] + std::sqrt(2.0 / referencePrecisionValues[j]) * nodeValues[k];
        double weightValue = weightValues[k];

        for (unsigned int l = 0;l < numSamples;++l)
        {
          if (!subsetValues[l])
            continue;

          mixtureModel = mixtureModels[l];
          meanValues = mixtureModel["mean"];
          precisionValues = mixtureModel["precision"];
          mixingValues = mixtureModel["mixing"];

          double workValue = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

          if (l == i)
            logValue -= (numSubset - 1.0) * workValue;
          else
            logValue += workValue;
        }

        logValue /= numSubset;

        logIntegral += weightValue * logValue;
        squareLogIntegral += weightValue * logValue * logValue;
      }

      totalLogIntegral += referenceMixingValue * logIntegral;
      totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
    }

    double workSquaredDistance = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;

    if (workSquaredDistance < referenceSquaredDistance || pos == 0)
    {
      referenceSquaredDistance = workSquaredDistance;
      referenceIndex = i;
    }

    ++pos;
  }

  // Then, get distances to mean wrt new reference
  referenceModel = mixtureModels[referenceIndex];
  referenceMeanValues = referenceModel["mean"];
  referencePrecisionValues = referenceModel["precision"];
  referenceMixingValues = referenceModel["mixing"];
  unsigned int referenceNumComponents = referenceModel.nrows();

  for (unsigned int i = 0;i < numSamples;++i)
  {
    if (i == referenceIndex)
    {
      squaredDistancesToMean[i] = referenceSquaredDistance;
      continue;
    }

    double totalLogIntegral = 0.0;
    double totalSquareLogIntegral = 0.0;

    for (unsigned int j = 0;j < referenceNumComponents;++j)
    {
      double logIntegral = 0.0;
      double squareLogIntegral = 0.0;
      double referenceMixingValue = referenceMixingValues[j];

      for (unsigned int k = 0;k < numNodes;++k)
      {
        double logValue = 0.0;
        double inputValue = referenceMeanValues[j] + std::sqrt(2.0 / referencePrecisionValues[j]) * nodeValues[k];
        double weightValue = weightValues[k];

        for (unsigned int l = 0;l < numSamples;++l)
        {
          if (!subsetValues[l])
            continue;

          mixtureModel = mixtureModels[l];
          meanValues = mixtureModel["mean"];
          precisionValues = mixtureModel["precision"];
          mixingValues = mixtureModel["mixing"];

          double workValue = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

          if (l == i)
            logValue += (numSubset - 1.0) * workValue;
          else
            logValue -= workValue;
        }

        logValue /= numSubset;

        if (!subsetValues[i])
        {
          mixtureModel = mixtureModels[i];
          meanValues = mixtureModel["mean"];
          precisionValues = mixtureModel["precision"];
          mixingValues = mixtureModel["mixing"];
          logValue += GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);
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
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numInputs = inputData.size();

  Rcpp::NumericVector outputValues(numInputs * (numInputs - 1) / 2);
  Rcpp::DataFrame firstModel, secondModel;
  Rcpp::NumericVector firstMeanValues, firstPrecisionValues, firstMixingValues;
  Rcpp::NumericVector secondMeanValues, secondPrecisionValues, secondMixingValues;
  std::vector <double> workVector, logValues(numInputs);

  // First define the reference Gaussian mixture wrt which integration will take place with GH.
  // The best for accuracy is to concatenate all input mixture models into one.
  // Note however that it is computationally expensive.

  std::vector<double> referenceMeanValues, referencePrecisionValues, referenceMixingValues;
  unsigned int numComponents = BuildReferenceModel(inputData, referenceMeanValues, referencePrecisionValues, referenceMixingValues);

  // Now integrate to get distances

  for (unsigned int i = 0;i < numInputs - 1;++i)
  {
    firstModel = inputData[i];
    firstMeanValues = firstModel["mean"];
    firstPrecisionValues = firstModel["precision"];
    firstMixingValues = firstModel["mixing"];

    for (unsigned int j = i + 1;j < numInputs;++j)
    {
      secondModel = inputData[j];
      secondMeanValues = secondModel["mean"];
      secondPrecisionValues = secondModel["precision"];
      secondMixingValues = secondModel["mixing"];

      double workScalar = GetSquaredDistance(
        firstMeanValues, firstPrecisionValues, firstMixingValues,
        secondMeanValues, secondPrecisionValues, secondMixingValues,
        referenceMeanValues, referencePrecisionValues, referenceMixingValues,
        nodeValues, weightValues, workVector);

      // workScalar += 1 * (firstModel.nrows() - secondModel.nrows()) * (firstModel.nrows() - secondModel.nrows());
      outputValues[numInputs * i - (i + 1) * i / 2 + j - i - 1] = workScalar;
    }
  }

  return outputValues;
}

double GetMeanNormalizationFactor(
    const Rcpp::List &inputData,
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

  for (unsigned int i = 0;i < numComponents;++i)
  {
    double mixingValue = referenceMixingValues[i];
    double integralValue = 0.0;

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double nodeValue = nodeValues[j];
      double weightValue = weightValues[j];
      double inputValue = referenceMeanValues[i] + std::sqrt(2.0 * referencePrecisionValues[i]) * nodeValue;

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

      integralValue += weightValue * std::exp(logarithmArithmeticMean - arithmeticMeanLogarithm);
    }

    totalIntegralValue += mixingValue * integralValue;
  }

  return totalIntegralValue;
}

// TO BE CHECKED
// ADD GetSquaredDistanceToMean()
double GetMeanExpectedValue(const Rcpp::List &inputData,
                            const Rcpp::NumericVector &nodeValues,
                            const Rcpp::NumericVector &weightValues)
{
  unsigned int numInputs = inputData.size();
  unsigned int numPoints = nodeValues.size();

  Rcpp::DataFrame inputModel;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector <double> workVector, logValues(numInputs);

  // First define the reference Gaussian mixture wrt which integration will take place with GH.
  // The best for accuracy is to concatenate all input mixture models into one.
  // Note however that it is computationally expensive.

  std::vector<double> referenceMeanValues, referencePrecisionValues, referenceMixingValues;
  unsigned int numComponents = BuildReferenceModel(inputData, referenceMeanValues, referencePrecisionValues, referenceMixingValues);

  // Now integrate

  double totalIntegralValueUp = 0.0;
  double totalIntegralValueDown = 0.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    double mixingValue = referenceMixingValues[i];
    double integralValueUp = 0.0;
    double integralValueDown = 0.0;

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double nodeValue = nodeValues[j];
      double weightValue = weightValues[j];
      double inputValue = referenceMeanValues[i] + std::sqrt(2.0 * referencePrecisionValues[i]) * nodeValue;
      double maxLogValue = 0.0;

      for (unsigned int k = 0;k < numInputs;++k)
      {
        inputModel = inputData[k];
        meanValues = inputModel["mean"];
        precisionValues = inputModel["precision"];
        mixingValues = inputModel["mixing"];

        double logValue = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, workVector);

        if (logValue > maxLogValue || k == 0)
          maxLogValue = logValue;

        logValues[k] = logValue;
      }

      double numeratorValue = 0.0;
      double denominatorValue = 0.0;

      for (unsigned int k = 0;k < numInputs;++k)
      {
        numeratorValue += logValues[k];
        denominatorValue += std::exp(logValues[k] - maxLogValue);
      }

      numeratorValue /= numInputs;
      numeratorValue -= maxLogValue;
      numeratorValue = std::exp(numeratorValue);

      integralValueUp += weightValue * inputValue * numeratorValue / denominatorValue;
      integralValueDown += weightValue * numeratorValue / denominatorValue;
    }

    totalIntegralValueUp += mixingValue * integralValueUp;
    totalIntegralValueDown += mixingValue * integralValueDown;
  }

  return totalIntegralValueUp / totalIntegralValueDown;
}

Rcpp::List PerformReassignment(
  const Rcpp::IntegerVector &inputMemberships,
  const Rcpp::List &inputModels,
  const Rcpp::IntegerVector &referenceModels,
  const Rcpp::NumericVector &nodeValues,
  const Rcpp::NumericVector &weightValues,
  const double alpha)
{
  // Mean is obtained from membership vector
  // Output is new membership vector for a specific group as defined by membershipVector
  unsigned int numSamples = inputModels.size();
  unsigned int numClusters = referenceModels.size();
  Rcpp::IntegerVector outputMemberships(numSamples);
  Rcpp::NumericVector outputDistances(numSamples);

  for (unsigned int i = 0;i < numSamples;++i)
  {
    unsigned int bestCluster = 0;
    double smallestDistance = 0.0;

    for (unsigned int j = 0;j < numClusters;++j)
    {
      double workDistance = GetDistanceToMean(inputModels, i, inputMemberships, j, referenceModels[j], nodeValues, weightValues, alpha);

      if (workDistance < smallestDistance || j == 0)
      {
        smallestDistance = workDistance;
        bestCluster = j;
      }
    }

    outputMemberships[i] = bestCluster;
    outputDistances[i] = smallestDistance;
  }

  return Rcpp::List::create(
    Rcpp::Named("memberships") = outputMemberships,
    Rcpp::Named("distances") = outputDistances
  );
}

Rcpp::IntegerVector ComputeNewReferences(
    const Rcpp::IntegerVector &inputReferences,
    const Rcpp::List &inputModels,
    const Rcpp::IntegerVector &membershipValues,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numSamples = inputModels.size();
  unsigned int numClusters = inputReferences.size();
  Rcpp::IntegerVector outputReferences(numClusters);

  for (unsigned int j = 0;j < numClusters;++j)
  {
    unsigned int bestReference = 0;
    unsigned int pos = 0;
    double smallestDistance = 0.0;

    for (unsigned int i = 0;i < numSamples;++i)
    {
      if (membershipValues[i] != j)
        continue;

      double workDistance = GetDistanceToMean(inputModels, i, membershipValues, j, inputReferences[j], nodeValues, weightValues);

      if (workDistance < smallestDistance || pos == 0)
      {
        smallestDistance = workDistance;
        bestReference = i;
      }

      ++pos;
    }

    outputReferences[j] = bestReference;
  }

  return outputReferences;
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

double GetDistanceToMean(
    const Rcpp::List &mixtureModels,
    const unsigned int observationIndex,
    const Rcpp::IntegerVector &membershipValues,
    const unsigned int clusterIndex,
    const unsigned int referenceIndex,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues,
    const double alpha)
{
  unsigned int numNodes = nodeValues.size();
  unsigned int numSamples = mixtureModels.size();

  unsigned int numSubset = 0;
  for (unsigned int i = 0;i < numSamples;++i)
  {
    if (membershipValues[i] == clusterIndex)
      ++numSubset;
  }

  unsigned int numUselessSamples = std::floor(alpha * numSubset);
  std::vector<double> sampleValues(numSubset);

  Rcpp::DataFrame mixtureModel;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector<double> logDensityValues;

  Rcpp::DataFrame referenceModel = Rcpp::as<Rcpp::DataFrame>(mixtureModels[referenceIndex]);
  Rcpp::NumericVector referenceMeanValues = referenceModel["mean"];
  Rcpp::NumericVector referencePrecisionValues = referenceModel["precision"];
  Rcpp::NumericVector referenceMixingValues = referenceModel["mixing"];
  unsigned int referenceNumComponents = referenceModel.nrows();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;

  for (unsigned int j = 0;j < referenceNumComponents;++j)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    double referenceMixingValue = referenceMixingValues[j];

    for (unsigned int k = 0;k < numNodes;++k)
    {
      double inputValue = referenceMeanValues[j] + std::sqrt(2.0 / referencePrecisionValues[j]) * nodeValues[k];
      double weightValue = weightValues[k];

      // Compute logValue of i-th observed mixture
      mixtureModel = mixtureModels[observationIndex];
      meanValues = mixtureModel["mean"];
      precisionValues = mixtureModel["precision"];
      mixingValues = mixtureModel["mixing"];

      double logValue1 = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

      // Compute logValue of alpha-trimmed mean
      unsigned int pos = 0;
      for (unsigned int l = 0;l < numSamples;++l)
      {
        if (membershipValues[l] != clusterIndex)
          continue;

        mixtureModel = mixtureModels[l];
        meanValues = mixtureModel["mean"];
        precisionValues = mixtureModel["precision"];
        mixingValues = mixtureModel["mixing"];

        sampleValues[pos] = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);
        ++pos;
      }

      std::nth_element(sampleValues.begin(), sampleValues.begin() + numUselessSamples, sampleValues.end());

      double logValue2 = 0.0;
      for (unsigned int l = numUselessSamples;l < numSubset;++l)
        logValue2 += sampleValues[l];

      logValue2 /= (numSubset - numUselessSamples);

      logIntegral += weightValue * (logValue1 - logValue2);
      squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
  }

  return totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
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
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues,
    std::vector<double> &workVector)
{
  unsigned int numPoints = nodeValues.size();
  unsigned int referenceNumComponents = referenceMeanValues.size();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;

  for (unsigned int j = 0;j < referenceNumComponents;++j)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    double referenceMixingValue = referenceMixingValues[j];

    for (unsigned int k = 0;k < numPoints;++k)
    {
      double inputValue = referenceMeanValues[j] + std::sqrt(2.0 / referencePrecisionValues[j]) * nodeValues[k];
      double weightValue = weightValues[k];

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
    for (unsigned int i = 0;i < internalNumComponents;++i)
      mixingValues[i] /= numInputs;
    referenceMixingValues.insert(referenceMixingValues.end(), mixingValues.begin(), mixingValues.end());

    numComponents += internalNumComponents;
  }

  return numComponents;
}
