#include "integrandClasses.h"

void GenericIntegrand::SetReferenceMixture(
    const VectorType &means,
    const VectorType &precisions,
    const VectorType &mixings,
    const unsigned int ncomponents)
{
  m_ReferenceMixture.SetMeanValues(means);
  m_ReferenceMixture.SetPrecisionValues(precisions);
  m_ReferenceMixture.SetMixingValues(mixings);
  m_ReferenceMixture.SetNumberOfComponents(ncomponents);
}

double LogDifference::operator() (const double& x) const
{
  double logValue1 = m_FirstMixture.GetLogDensity(x);
  double logValue2 = m_SecondMixture.GetLogDensity(x);
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  return (logValue1 - logValue2) * std::exp(logValueRef);
}

double SquaredLogDifference::operator() (const double& x) const
{
  double logValue1 = m_FirstMixture.GetLogDensity(x);
  double logValue2 = m_SecondMixture.GetLogDensity(x);
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  return (logValue1 - logValue2) * (logValue1 - logValue2) * std::exp(logValueRef);
}

double LogDifferenceFirstDerivative::operator() (const double& x) const
{
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  double firstDerivative = m_FirstMixture.GetLogDensityShiftDerivative(x);
  return firstDerivative * std::exp(logValueRef);
}

double SquaredLogDifferenceFirstDerivative::operator() (const double& x) const
{
  double logValue1 = m_FirstMixture.GetLogDensity(x);
  double logValue2 = m_SecondMixture.GetLogDensity(x);
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  double firstDerivative = m_FirstMixture.GetLogDensityShiftDerivative(x);
  return firstDerivative * (logValue1 - logValue2) * std::exp(logValueRef);
}

double LogDifferenceSecondDerivative::operator() (const double& x) const
{
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  double secondDerivative = -m_SecondMixture.GetLogDensityShiftDerivative(x);
  return secondDerivative * std::exp(logValueRef);
}

double SquaredLogDifferenceSecondDerivative::operator() (const double& x) const
{
  double logValue1 = m_FirstMixture.GetLogDensity(x);
  double logValue2 = m_SecondMixture.GetLogDensity(x);
  double logValueRef = m_ReferenceMixture.GetLogDensity(x);
  double secondDerivative = -m_SecondMixture.GetLogDensityShiftDerivative(x);
  return secondDerivative * (logValue1 - logValue2) * std::exp(logValueRef);
}
