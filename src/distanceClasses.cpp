#include "distanceClasses.h"
#include <cmath>

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
  // First, reinitialize reference
  m_ReferenceMeanValues.clear();
  m_ReferencePrecisionValues.clear();
  m_ReferenceMixingValues.clear();
  m_ReferenceNumberOfComponents = 0;

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
  unsigned int firstNumComponents = m_FirstMixture.GetNumberOfComponents();
  unsigned int secondNumComponents = m_SecondMixture.GetNumberOfComponents();
  m_ReferenceNumberOfComponents = firstNumComponents + secondNumComponents;
  m_ReferenceMeanValues.resize(m_ReferenceNumberOfComponents);
  m_ReferencePrecisionValues.resize(m_ReferenceNumberOfComponents);
  m_ReferenceMixingValues.resize(m_ReferenceNumberOfComponents);

  m_WorkVector = m_SecondMixture.GetMeanValues();
  m_ReferenceMeanValues.insert(m_ReferenceMeanValues.begin() + firstNumComponents, m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_SecondMixture.GetPrecisionValues();
  m_ReferencePrecisionValues.insert(m_ReferencePrecisionValues.begin() + firstNumComponents, m_WorkVector.begin(), m_WorkVector.end());

  m_WorkVector = m_SecondMixture .GetMixingValues();
  for (unsigned int i = 0;i < secondNumComponents;++i)
    m_WorkVector[i] /= 2.0;
  m_ReferenceMixingValues.insert(m_ReferenceMixingValues.begin() + firstNumComponents, m_WorkVector.begin(), m_WorkVector.end());
}

double SquaredBayesDistance::EvaluateSquaredDistance()
{
  unsigned int numPoints = m_QuadraturePoints.size();

  double totalSLDIntegral = 0.0;
  double totalLDIntegral = 0.0;
  double totalJLDIntegral = 0.0;
  double totalSLDJIntegral = 0.0;
  double totalLDJIntegral = 0.0;
  double totalJIntegral = 0.0;

  for (unsigned int i = 0;i < m_ReferenceNumberOfComponents;++i)
  {
    double refMeanValue = (i < m_FirstMixture.GetNumberOfComponents()) ? m_ReferenceMeanValues[i] : m_SecondMixtureCopy.GetMeanValues()[i-m_FirstMixture.GetNumberOfComponents()];

    double  sldIntegral = 0.0;
    double   ldIntegral = 0.0;
    double  jldIntegral = 0.0;
    double sldjIntegral = 0.0;
    double  ldjIntegral = 0.0;
    double    jIntegral = 0.0;

    double referenceWeight = m_ReferenceMixingValues[i];

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double quadWeight = m_QuadratureWeights[j];
      double quadPoint = refMeanValue + std::sqrt(2.0 / m_ReferencePrecisionValues[i]) * m_QuadraturePoints[j];

      double logValue1 = m_FirstMixture.GetLogDensity(quadPoint);
      double logValue2 = m_SecondMixtureCopy.GetLogDensity(quadPoint);
      double logDifference = logValue1 - logValue2;
      double shiftLDJacobian = m_SecondMixtureCopy.GetLogDensityShiftDerivative(quadPoint);

      sldIntegral += quadWeight * logDifference * logDifference;
      ldIntegral += quadWeight * logDifference;
      jldIntegral += quadWeight * logDifference * shiftLDJacobian;
      if (i >= m_FirstMixture.GetNumberOfComponents())
      {
        double workScalar = quadWeight * logDifference * (quadPoint - refMeanValue);
        sldjIntegral += logDifference * workScalar;
        ldjIntegral += workScalar;
      }
      jIntegral += quadWeight * shiftLDJacobian;
    }

    totalSLDIntegral += referenceWeight * sldIntegral;
    totalLDIntegral += referenceWeight * ldIntegral;
    totalJLDIntegral += referenceWeight * jldIntegral;
    if (i >= m_FirstMixture.GetNumberOfComponents())
    {
      totalSLDJIntegral += m_ReferencePrecisionValues[i] * referenceWeight * sldjIntegral;
      totalLDJIntegral += m_ReferencePrecisionValues[i] * referenceWeight * ldjIntegral;
    }
    totalJIntegral += referenceWeight * jIntegral;
  }

  double costValue = totalSLDIntegral / std::sqrt(M_PI) - totalLDIntegral * totalLDIntegral / M_PI;

  m_SecondShiftDerivative = -2.0 * totalJLDIntegral;
  m_SecondShiftDerivative += totalSLDJIntegral;
  m_SecondShiftDerivative /= std::sqrt(M_PI);
  m_SecondShiftDerivative -= 2.0 * totalLDIntegral * (totalLDJIntegral - totalJIntegral) / M_PI;

  return costValue;
}
