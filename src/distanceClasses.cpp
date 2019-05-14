#include "distanceClasses.h"
#include <cmath>

double SquaredBayesDistance::f_grad(Numer::Constvec& x, Numer::Refvec grad)
{
  m_FirstMixtureCopy = m_FirstMixture;
  m_WorkVector = m_FirstMixture.GetMeanValues();
  for (unsigned int i = 0;i < m_FirstMixture.GetNumberOfComponents();++i)
    m_WorkVector[i] += x[0];
  m_FirstMixtureCopy.SetMeanValues(m_WorkVector);

  m_SecondMixtureCopy = m_SecondMixture;
  m_WorkVector = m_SecondMixture.GetMeanValues();
  for (unsigned int i = 0;i < m_SecondMixture.GetNumberOfComponents();++i)
    m_WorkVector[i] += x[1];
  m_SecondMixtureCopy.SetMeanValues(m_WorkVector);

  double costValue = this->EvaluateSquaredDistance();

  grad[0] = m_FirstShiftDerivative;
  grad[1] = m_SecondShiftDerivative;

  if (!std::isfinite(costValue))
    Rcpp::Rcout << costValue << std::endl;
  if (!std::isfinite(grad[0]) || !std::isfinite(grad[1]))
    Rcpp::Rcout << grad << std::endl;

  return costValue;
}

void SquaredBayesDistance::SetDataPoints(const Rcpp::List &inputData)
{
  m_TotalNumberOfComponents = this->SetReferenceModel(inputData);
  m_DataPoints = inputData;
}

void SquaredBayesDistance::SetInputValues(const unsigned int input, const unsigned int index)
{
  m_WorkTibble = m_DataPoints[index];

  if (input == 0)
  {
    m_FirstMixture.SetInput(m_WorkTibble);
    return;
  }
  m_SecondMixture.SetInput(m_WorkTibble);
}

unsigned int SquaredBayesDistance::SetReferenceModel(const Rcpp::List &inputData)
{
  unsigned int numInputs = inputData.size();

  m_ReferenceMeanValues.clear();
  m_ReferencePrecisionValues.clear();
  m_ReferenceMixingValues.clear();

  m_ReferenceMeanValue = 0.0;
  unsigned int numComponents = 0;

  for (unsigned int i = 0;i < numInputs;++i)
  {
    m_WorkTibble = inputData[i];
    m_WorkMixture.SetInput(m_WorkTibble);

    unsigned int internalNumComponents = m_WorkMixture.GetNumberOfComponents();

    m_WorkVector = m_WorkMixture.GetMeanValues();
    m_ReferenceMeanValues.insert(m_ReferenceMeanValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    m_WorkVector = m_WorkMixture.GetPrecisionValues();
    m_ReferencePrecisionValues.insert(m_ReferencePrecisionValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    m_WorkVector = m_WorkMixture.GetMixingValues();
    for (unsigned int j = 0;j < internalNumComponents;++j)
      m_WorkVector[j] /= numInputs;
    m_ReferenceMixingValues.insert(m_ReferenceMixingValues.end(), m_WorkVector.begin(), m_WorkVector.end());

    m_ReferenceMeanValue += m_WorkMixture.GetMean();
    numComponents += internalNumComponents;
  }

  m_ReferenceMeanValue /= numInputs;

  return numComponents;
}

double SquaredBayesDistance::EvaluateSquaredDistance()
{
  unsigned int numPoints = m_QuadraturePoints.size();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;
  // Other four integrals for computing Jacobian
  double totalLogIntegralWRTFirst = 0.0;
  double totalSquareLogIntegralWRTFirst = 0.0;
  double totalLogIntegralWRTSecond = 0.0;
  double totalSquareLogIntegralWRTSecond = 0.0;

  for (unsigned int i = 0;i < m_TotalNumberOfComponents;++i)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    // Other four integrals for computing Jacobian
    double logIntegralWRTFirst = 0.0;
    double squareLogIntegralWRTFirst = 0.0;
    double logIntegralWRTSecond = 0.0;
    double squareLogIntegralWRTSecond = 0.0;

    double referenceMixingValue = m_ReferenceMixingValues[i];

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double inputValue = m_ReferenceMeanValues[i] + std::sqrt(2.0 / m_ReferencePrecisionValues[i]) * m_QuadraturePoints[j];
      double weightValue = m_QuadratureWeights[j];

      double logValue1 = m_FirstMixtureCopy.GetLogDensity(inputValue);
      double logValue2 = m_SecondMixtureCopy.GetLogDensity(inputValue);

      logIntegral += weightValue * (logValue1 - logValue2);
      squareLogIntegral += weightValue * (logValue1 - logValue2) * (logValue1 - logValue2);

      double firstJacobianValue = m_FirstMixtureCopy.GetLogDensityShiftDerivative(inputValue);
      logIntegralWRTFirst += weightValue * firstJacobianValue;
      squareLogIntegralWRTFirst += weightValue * firstJacobianValue * (logValue1 - logValue2);

      double secondJacobianValue = -m_SecondMixtureCopy.GetLogDensityShiftDerivative(inputValue);
      logIntegralWRTSecond += weightValue * secondJacobianValue;
      squareLogIntegralWRTSecond += weightValue * secondJacobianValue * (logValue1 - logValue2);
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;

    totalLogIntegralWRTFirst += referenceMixingValue * logIntegralWRTFirst;
    totalSquareLogIntegralWRTFirst += referenceMixingValue * squareLogIntegralWRTFirst;

    totalLogIntegralWRTSecond += referenceMixingValue * logIntegralWRTSecond;
    totalSquareLogIntegralWRTSecond += referenceMixingValue * squareLogIntegralWRTSecond;
  }

  double costValue = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;

  m_FirstShiftDerivative = totalSquareLogIntegralWRTFirst / std::sqrt(M_PI);
  m_FirstShiftDerivative -= totalLogIntegralWRTFirst * totalLogIntegral / M_PI;
  m_FirstShiftDerivative *= 2.0;

  m_SecondShiftDerivative = totalSquareLogIntegralWRTSecond / std::sqrt(M_PI);
  m_SecondShiftDerivative -= totalLogIntegralWRTSecond * totalLogIntegral / M_PI;
  m_SecondShiftDerivative *= 2.0;

  return costValue;
}
