#include "distanceClasses.h"
#include <cmath>

void SquaredBayesDistance::InitializeReferenceModel()
{
  m_ReferenceMeanValues.clear();
  m_ReferencePrecisionValues.clear();
  m_ReferenceMixingValues.clear();
  m_ReferenceNumberOfComponents = 0;
}

double SquaredBayesDistance::f_grad(Numer::Constvec& x, Numer::Refvec grad)
{
  m_SecondMixtureCopy = m_SecondMixture;
  m_WorkVector = m_SecondMixture.GetMeanValues();
  for (unsigned int i = 0;i < m_SecondMixture.GetNumberOfComponents();++i)
    m_WorkVector[i] += x[0];
  m_SecondMixtureCopy.SetMeanValues(m_WorkVector);

  double costValue = this->EvaluateSquaredDistance();

  grad[0] = m_SecondShiftDerivative;

  if (!std::isfinite(costValue))
    Rcpp::Rcout << costValue << std::endl;
  if (!std::isfinite(grad[0]))
    Rcpp::Rcout << grad << std::endl;

  return costValue;
}

void SquaredBayesDistance::SetInput1(const Rcpp::DataFrame &x)
{
  m_FirstMixture.SetInput(x);

  unsigned int numComponents = m_FirstMixture.GetNumberOfComponents();
  m_ReferenceNumberOfComponents += numComponents;

  m_WorkVector = m_FirstMixture.GetMeanValues();
  m_ReferenceMeanValues.insert(m_ReferenceMeanValues.end(), m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_FirstMixture.GetPrecisionValues();
  m_ReferencePrecisionValues.insert(m_ReferencePrecisionValues.end(), m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_FirstMixture.GetMixingValues();
  for (unsigned int i = 0;i < numComponents;++i)
    m_WorkVector[i] /= 2.0;
  m_ReferenceMixingValues.insert(m_ReferenceMixingValues.end(), m_WorkVector.begin(), m_WorkVector.end());
}

void SquaredBayesDistance::SetInput2(const Rcpp::DataFrame &x)
{
  m_SecondMixture.SetInput(x);

  unsigned int numComponents = m_SecondMixture.GetNumberOfComponents();
  m_ReferenceNumberOfComponents += numComponents;

  m_WorkVector = m_SecondMixture.GetMeanValues();
  m_ReferenceMeanValues.insert(m_ReferenceMeanValues.end(), m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_SecondMixture.GetPrecisionValues();
  m_ReferencePrecisionValues.insert(m_ReferencePrecisionValues.end(), m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_SecondMixture .GetMixingValues();
  for (unsigned int i = 0;i < numComponents;++i)
    m_WorkVector[i] /= 2.0;
  m_ReferenceMixingValues.insert(m_ReferenceMixingValues.end(), m_WorkVector.begin(), m_WorkVector.end());
}

double SquaredBayesDistance::EvaluateSquaredDistance()
{
  unsigned int numPoints = m_QuadraturePoints.size();

  double totalLogIntegral = 0.0;
  double totalSquareLogIntegral = 0.0;
  // Other two integrals for computing Jacobian
  double totalLogIntegralWRTSecond = 0.0;
  double totalSquareLogIntegralWRTSecond = 0.0;

  for (unsigned int i = 0;i < m_ReferenceNumberOfComponents;++i)
  {
    double logIntegral = 0.0;
    double squareLogIntegral = 0.0;
    // Other two integrals for computing Jacobian
    double logIntegralWRTSecond = 0.0;
    double squareLogIntegralWRTSecond = 0.0;

    double referenceMixingValue = m_ReferenceMixingValues[i];

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double weightValue = m_QuadratureWeights[j];
      double inputValue = m_ReferenceMeanValues[i] + std::sqrt(2.0 / m_ReferencePrecisionValues[i]) * m_QuadraturePoints[j];

      double logValue1 = m_FirstMixture.GetLogDensity(inputValue);
      double logValue2 = m_SecondMixtureCopy.GetLogDensity(inputValue);
      double logDifference = logValue1 - logValue2;

      logIntegral += weightValue * logDifference;
      squareLogIntegral += weightValue * logDifference * logDifference;

      double secondJacobianValue = -m_SecondMixtureCopy.GetLogDensityShiftDerivative(inputValue);
      logIntegralWRTSecond += weightValue * secondJacobianValue;
      squareLogIntegralWRTSecond += weightValue * secondJacobianValue * logDifference;
    }

    totalLogIntegral += referenceMixingValue * logIntegral;
    totalSquareLogIntegral += referenceMixingValue * squareLogIntegral;

    totalLogIntegralWRTSecond += referenceMixingValue * logIntegralWRTSecond;
    totalSquareLogIntegralWRTSecond += referenceMixingValue * squareLogIntegralWRTSecond;
  }

  double costValue = totalSquareLogIntegral / std::sqrt(M_PI) - totalLogIntegral * totalLogIntegral / M_PI;

  m_SecondShiftDerivative = totalSquareLogIntegralWRTSecond / std::sqrt(M_PI);
  m_SecondShiftDerivative -= totalLogIntegralWRTSecond * totalLogIntegral / M_PI;
  m_SecondShiftDerivative *= 2.0;

  return costValue;
}
