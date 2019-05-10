#include "distanceClasses.h"

double SquaredBayesDistance::f_grad(Constvec& x, Refvec grad)
{
  unsigned firstMeanSize = m_FirstMeanValues.size();
  VectorType firstMeanValues(firstMeanSize);
  for (unsigned int i = 0;i < firstMeanSize;++i)
    firstMeanValues[i] = m_FirstMeanValues[i] + x[0];

  VectorType referenceMeanValues(m_TotalNumberOfComponents);
  for (unsigned int i = 0;i < m_TotalNumberOfComponents;++i)
    referenceMeanValues[i] = m_ReferenceMeanValues[i] + x[1];

  double jacobianValueWRTFirst, jacobianValueWRTReference;

  double costValue = this->EvaluateSquaredDistance(
    firstMeanValues, m_FirstPrecisionValues, m_FirstMixingValues,
    m_SecondMeanValues, m_SecondPrecisionValues, m_SecondMixingValues,
    referenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues,
    jacobianValueWRTFirst, jacobianValueWRTReference);

  grad[0] = jacobianValueWRTFirst;
  grad[1] = jacobianValueWRTReference;

  // Rcpp::Rcout << costValue << std::endl;
  // Rcpp::Rcout << grad << std::endl;

  return costValue;
}

void SquaredBayesDistance::SetDataPoints(const Rcpp::List &inputData)
{
  // Build reference model
  m_TotalNumberOfComponents = this->SetReferenceModel(inputData);
  m_DataPoints = inputData;
}

void SquaredBayesDistance::SetInputValues(const unsigned int input, const unsigned int index)
{
  m_WorkTibble = m_DataPoints[index];
  if (input == 0)
  {
    m_FirstMeanValues = Rcpp::as<VectorType>(m_WorkTibble["mean"]);
    m_FirstPrecisionValues = Rcpp::as<VectorType>(m_WorkTibble["precision"]);
    m_FirstMixingValues = Rcpp::as<VectorType>(m_WorkTibble["mixing"]);
    return;
  }
  m_SecondMeanValues = Rcpp::as<VectorType>(m_WorkTibble["mean"]);
  m_SecondPrecisionValues = Rcpp::as<VectorType>(m_WorkTibble["precision"]);
  m_SecondMixingValues = Rcpp::as<VectorType>(m_WorkTibble["mixing"]);
}

unsigned int SquaredBayesDistance::SetReferenceModel(const Rcpp::List &inputData)
{
  unsigned int numInputs = inputData.size();

  m_ReferenceMeanValues.clear();
  m_ReferencePrecisionValues.clear();
  m_ReferenceMixingValues.clear();

  unsigned int numComponents = 0;

  for (unsigned int i = 0;i < numInputs;++i)
  {
    m_WorkTibble = inputData[i];
    unsigned int internalNumComponents = m_WorkTibble.nrows();

    m_WorkVector = Rcpp::as<VectorType>(m_WorkTibble["mean"]);
    m_ReferenceMeanValues.insert(m_ReferenceMeanValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    m_WorkVector = Rcpp::as<VectorType>(m_WorkTibble["precision"]);
    m_ReferencePrecisionValues.insert(m_ReferencePrecisionValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    m_WorkVector = Rcpp::as<VectorType>(m_WorkTibble["mixing"]);
    for (unsigned int j = 0;j < internalNumComponents;++j)
      m_WorkVector[j] /= numInputs;
    m_ReferenceMixingValues.insert(m_ReferenceMixingValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    numComponents += internalNumComponents;
  }

  return numComponents;
}

double SquaredBayesDistance::EvaluateSquaredDistance(
    const VectorType &firstMeanValues, const VectorType &firstPrecisionValues, const VectorType &firstMixingValues,
    const VectorType &secondMeanValues, const VectorType &secondPrecisionValues, const VectorType &secondMixingValues,
    const VectorType &referenceMeanValues, const VectorType &referencePrecisionValues, const VectorType &referenceMixingValues,
    double &jacobianValueWRTFirst, double &jacobianValueWRTReference)
{
  unsigned int numPoints = m_QuadraturePoints.size();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;
  // Other four integrals for computing Jacobian
  double totalLogIntegralWRTFirst = 0.0;
  double totalSquareLogIntegralWRTFirst = 0.0;
  double totalLogIntegralWRTReference = 0.0;
  double totalSquareLogIntegralWRTReference = 0.0;

  for (unsigned int i = 0;i < m_TotalNumberOfComponents;++i)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    // Other four integrals for computing Jacobian
    double logIntegralWRTFirst = 0.0;
    double squareLogIntegralWRTFirst = 0.0;
    double logIntegralWRTReference = 0.0;
    double squareLogIntegralWRTReference = 0.0;
    double referenceMixingValue = referenceMixingValues[i];
    // Other weight for some integrals involved in Jacobian
    double referenceMixingValueWRTReference = referenceMixingValue * referencePrecisionValues[i];

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double inputValue = referenceMeanValues[i] + std::sqrt(2.0 / referencePrecisionValues[i]) * m_QuadraturePoints[j];
      double weightValue = m_QuadratureWeights[j];

      double logValue1 = this->EvaluateLogDensity(inputValue, firstMeanValues, firstPrecisionValues, firstMixingValues);
      double logValue2 = this->EvaluateLogDensity(inputValue, secondMeanValues, secondPrecisionValues, secondMixingValues);

      logIntegral += weightValue * (logValue1 - logValue2);
      squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);

      double firstJacobianValue = this->EvaluateLogDensityJacobian(inputValue, firstMeanValues, firstPrecisionValues, firstMixingValues);
      logIntegralWRTFirst += weightValue * firstJacobianValue;
      squareLogIntegralWRTFirst += weightValue * firstJacobianValue * (logValue1 - logValue2);

      logIntegralWRTReference += weightValue * (logValue1 - logValue2) * (inputValue - referenceMeanValues[i]);
      squareLogIntegralWRTReference += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2) * (inputValue - referenceMeanValues[i]);
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;

    totalLogIntegralWRTFirst += referenceMixingValue * logIntegralWRTFirst;
    totalSquareLogIntegralWRTFirst += referenceMixingValue * squareLogIntegralWRTFirst;
    totalLogIntegralWRTReference += referenceMixingValueWRTReference * logIntegralWRTReference;
    totalSquareLogIntegralWRTReference += referenceMixingValueWRTReference * squareLogIntegralWRTReference;
  }

  double costValue = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;

  jacobianValueWRTFirst = totalSquareLogIntegralWRTFirst / std::sqrt(M_PI);
  jacobianValueWRTFirst -= totalLogIntegralWRTFirst * totalLogIntegral / M_PI;
  jacobianValueWRTFirst *= 2.0;

  jacobianValueWRTReference = totalSquareLogIntegralWRTReference / std::sqrt(M_PI);
  jacobianValueWRTReference -= 2.0 * totalLogIntegralWRTReference * totalLogIntegral / M_PI;

  return costValue;
}

double SquaredBayesDistance::EvaluateLogDensity(
    const double inputValue,
    const VectorType &meanValues,
    const VectorType &precisionValues,
    const VectorType &mixingValues)
{
  unsigned int numComponents = meanValues.size();
  m_WorkVector.resize(numComponents);

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

    m_WorkVector[i] = 0.5 * std::log(precisionValues[i] / (2.0 * M_PI)) - workValue / 2.0;
  }

  double insideLogValue = 1.0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    if (i == indexOfMinimalExponent)
      continue;

    insideLogValue += std::exp(m_WorkVector[i] - m_WorkVector[indexOfMinimalExponent]) * mixingValues[i] / mixingValues[indexOfMinimalExponent];
  }

  // Compute log of mixture density
  return std::log(mixingValues[indexOfMinimalExponent]) + m_WorkVector[indexOfMinimalExponent] + std::log(insideLogValue);
}

double SquaredBayesDistance::EvaluateLogDensityJacobian(
    const double inputValue,
    const VectorType &meanValues,
    const VectorType &precisionValues,
    const VectorType &mixingValues)
{
  unsigned int numComponents = meanValues.size();
  m_WorkVector.resize(numComponents);

  double maximalExponent = 0.0;
  unsigned int indexOfMaximalExponent = 0;

  for (unsigned int i = 0;i < numComponents;++i)
  {
    double workScalar = -precisionValues[i] * (inputValue - meanValues[i]) * (inputValue - meanValues[i]) / 2.0;

    if (workScalar > maximalExponent || i == 0)
    {
      maximalExponent = workScalar;
      indexOfMaximalExponent = i;
    }

    m_WorkVector[i] = workScalar;
  }

  double sumWeights = 0.0;
  for (unsigned int i = 0;i < numComponents;++i)
  {
    double workScalar = mixingValues[i] * std::sqrt(precisionValues[i] / (2.0 * M_PI)) * std::exp(m_WorkVector[i] - maximalExponent);
    m_WorkVector[i] = workScalar;
    sumWeights += workScalar;
  }

  if (sumWeights < 1.0e-16)
  {
    sumWeights = (double)numComponents;
    for (unsigned int i = 0;i < numComponents;++i)
      m_WorkVector[i] = 1.0;
  }

  double outputValue = 0.0;

  for (unsigned int i = 0;i < numComponents;++i)
    outputValue += precisionValues[i] * (inputValue - meanValues[i]) * m_WorkVector[i];

  outputValue /= sumWeights;

  return outputValue;
}
