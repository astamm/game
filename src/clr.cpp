#include "clr.h"

Rcpp::NumericVector GetLogDensityWRTLebesgue(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues)
{
  unsigned int numPoints = inputValues.size();
  Rcpp::NumericVector outputValues(numPoints);
  std::vector<double> logDensityValues;

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double inputValue = inputValues[i];
    outputValues[i] = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);
  }

  return outputValues;
}

Rcpp::NumericVector GetLogDensityWRTGaussian(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    const double referenceMean,
    const double referencePrecision)
{
  Rcpp::NumericVector lebesgueValues = GetLogDensityWRTLebesgue(inputValues, meanValues, precisionValues, mixingValues);
  Rcpp::NumericVector outputValues(inputValues.size());

  for (unsigned int i = 0;i < inputValues.size();++i)
  {
    double inputValue = inputValues[i];
    double workValue = referencePrecision * (inputValue - referenceMean) * (inputValue - referenceMean);
    double referenceValue = 0.5 * std::log(referencePrecision / (2.0 * M_PI)) - workValue / 2.0;
    outputValues[i] = lebesgueValues[i] - referenceValue;
  }

  return outputValues;
}

Rcpp::NumericVector GetCenteredLogRatio(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    const double referenceMean,
    const double referencePrecision,
    const double centeringValue)
{
  Rcpp::NumericVector gaussianValues = GetLogDensityWRTGaussian(inputValues, meanValues, precisionValues, mixingValues, referenceMean, referencePrecision);
  Rcpp::NumericVector outputValues(inputValues.size());

  for (unsigned int i = 0;i < inputValues.size();++i)
    outputValues[i] = gaussianValues[i] - centeringValue;

  return outputValues;
}

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

Rcpp::NumericVector GetMeanSquaredNorms(
    const Rcpp::List &mixtureModels,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numNodes = nodeValues.size();
  unsigned int numSamples = mixtureModels.size();
  Rcpp::NumericVector squaredNorms(numSamples);
  Rcpp::DataFrame referenceMixtureModel, mixtureModel;
  Rcpp::NumericVector referenceMeanValues, referencePrecisionValues, referenceMixingValues;
  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector<double> logDensityValues;

  for (unsigned int i = 0;i < numSamples;++i)
  {
    referenceMixtureModel = mixtureModels[i];
    referenceMeanValues = referenceMixtureModel["mean"];
    referencePrecisionValues = referenceMixtureModel["precision"];
    referenceMixingValues = referenceMixtureModel["mixing"];
    unsigned int referenceNumComponents = referenceMixtureModel.nrows();

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
          mixtureModel = mixtureModels[l];
          meanValues = mixtureModel["mean"];
          precisionValues = mixtureModel["precision"];
          mixingValues = mixtureModel["mixing"];

          double workValue = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

          if (l == i)
            logValue -= (numSamples - 1.0) * workValue;
          else
            logValue += workValue;
        }

        logValue /= numSamples;

        logIntegral += weightValue * logValue;
        squareLogIntegral += weightValue * logValue * logValue;
      }

      totalLogIntegral += referenceMixingValue * logIntegral;
      totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
    }

    squaredNorms[i] = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
  }

  return squaredNorms;
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
    const Rcpp::List &inputModels,
    const Rcpp::DataFrame &referenceModel,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numSamples = inputModels.size();
  Rcpp::NumericVector squaredDistances(numSamples * (numSamples - 1) / 2);

  Rcpp::NumericVector referenceMeanValues = referenceModel["mean"];
  Rcpp::NumericVector referencePrecisionValues = referenceModel["precision"];
  Rcpp::NumericVector referenceMixingValues = referenceModel["mixing"];
  unsigned int referenceNumComponents = referenceModel.nrows();

  unsigned int numNodes = nodeValues.size();
  std::vector<double> logDensityValues;

  Rcpp::DataFrame firstModel, secondModel;
  Rcpp::NumericVector firstMeanValues, firstPrecisionValues, firstMixingValues;
  Rcpp::NumericVector secondMeanValues, secondPrecisionValues, secondMixingValues;

  for (unsigned int i = 0;i < numSamples - 1;++i)
  {
    firstModel = inputModels[i];
    firstMeanValues = firstModel["mean"];
    firstPrecisionValues = firstModel["precision"];
    firstMixingValues = firstModel["mixing"];

    for (unsigned int j = i + 1;j < numSamples;++j)
    {
      secondModel = inputModels[j];
      secondMeanValues = secondModel["mean"];
      secondPrecisionValues = secondModel["precision"];
      secondMixingValues = secondModel["mixing"];

      double totalLogIntegral = 0.0;
      double totalSquareLogIntegral = 0.0;

      for (unsigned int k = 0;k < referenceNumComponents;++k)
      {
        double logIntegral = 0.0;
        double squareLogIntegral = 0.0;
        double referenceMixingValue = referenceMixingValues[k];

        for (unsigned int l = 0;l < numNodes;++l)
        {
          double inputValue = referenceMeanValues[k] + std::sqrt(2.0 / referencePrecisionValues[k]) * nodeValues[l];
          double weightValue = weightValues[l];

          double logValue1 = GetLogDensityValue(inputValue, firstMeanValues, firstPrecisionValues, firstMixingValues, logDensityValues);
          double logValue2 = GetLogDensityValue(inputValue, secondMeanValues, secondPrecisionValues, secondMixingValues, logDensityValues);

          logIntegral += weightValue * (logValue1 - logValue2);
          squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);
        }

        totalLogIntegral += referenceMixingValue * logIntegral;
        totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
      }

      double sqDistance = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
      // sqDistance += (firstModel.nrows() - secondModel.nrows()) * (firstModel.nrows() - secondModel.nrows());
      squaredDistances[numSamples * i - (i + 1) * i / 2 + j - i - 1] = sqDistance;
    }
  }

  return squaredDistances;
}

double GetLogDensityValue(
    const double inputValue,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    std::vector<double> &logDensityValues)
{
  unsigned int numComponents = meanValues.size();
  logDensityValues.resize(numComponents);

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

    logDensityValues[i] = 0.5 * std::log(precisionValues[i] / (2.0 * M_PI)) - workValue / 2.0;
  }

  double insideLogValue = 1.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    if (i == indexOfMinimalExponent)
      continue;

    insideLogValue += std::exp(logDensityValues[i] - logDensityValues[indexOfMinimalExponent]) * mixingValues[i] / mixingValues[indexOfMinimalExponent];
  }

  // Compute log of mixture density
  return std::log(mixingValues[indexOfMinimalExponent]) + logDensityValues[indexOfMinimalExponent] + std::log(insideLogValue);
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

double GetDistance(
    const Rcpp::DataFrame &firstModel,
    const Rcpp::DataFrame &secondModel,
    const Rcpp::DataFrame &referenceModel,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues)
{
  unsigned int numNodes = nodeValues.size();

  Rcpp::NumericVector meanValues, precisionValues, mixingValues;
  std::vector<double> logDensityValues;

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

      meanValues = firstModel["mean"];
      precisionValues = firstModel["precision"];
      mixingValues = firstModel["mixing"];

      double logValue1 = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

      meanValues = secondModel["mean"];
      precisionValues = secondModel["precision"];
      mixingValues = secondModel["mixing"];

      double logValue2 = GetLogDensityValue(inputValue, meanValues, precisionValues, mixingValues, logDensityValues);

      logIntegral += weightValue * (logValue1 - logValue2);
      squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;
  }

  return totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;
}


