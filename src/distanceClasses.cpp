#include "distanceClasses.h"
#include <cmath>

double SquaredBayesDistance::operator()(const ParametersType &x, ParametersType &grad)
{
    m_WorkMixture = m_SecondMixture;
    m_WorkMeanValues = m_SecondMeanValues;
    for (unsigned int i = 0;i < m_SecondNumberOfComponents;++i)
        m_WorkMeanValues[i] += x[0];
    m_WorkMixture.SetMeanValues(m_WorkMeanValues);

    double costValue = this->ComputeSquaredDistance();

    grad[0] = m_ShiftDerivative;

    if (!std::isfinite(costValue))
        Rcpp::Rcout << costValue << std::endl;
    if (!std::isfinite(grad[0]))
        Rcpp::Rcout << grad << std::endl;

    return costValue;
}

void SquaredBayesDistance::SetInput1(const Rcpp::DataFrame &x)
{
  m_FirstMixture.SetInput(x);
  m_FirstNumberOfComponents = m_FirstMixture.GetNumberOfComponents();
  m_FirstMeanValues = m_FirstMixture.GetMeanValues();
  m_FirstPrecisionValues = m_FirstMixture.GetPrecisionValues();
  m_FirstMixingValues = m_FirstMixture.GetMixingValues();
}

void SquaredBayesDistance::SetInput2(const Rcpp::DataFrame &x)
{
  m_SecondMixture.SetInput(x);
  m_SecondNumberOfComponents = m_SecondMixture.GetNumberOfComponents();
  m_SecondMeanValues = m_SecondMixture.GetMeanValues();
  m_SecondPrecisionValues = m_SecondMixture.GetPrecisionValues();
  m_SecondMixingValues = m_SecondMixture.GetMixingValues();
}

double SquaredBayesDistance::ComputeSquaredDistance()
{
  unsigned int numPoints = m_QuadraturePoints.size();

  double totalSLDIntegral = 0.0;
  double totalLDIntegral = 0.0;
  double totalJLDIntegral = 0.0;
  double totalSLDJIntegral = 0.0;
  double totalLDJIntegral = 0.0;
  double totalJIntegral = 0.0;

  for (unsigned int i = 0;i < m_FirstNumberOfComponents + m_SecondNumberOfComponents;++i)
  {
    double referenceMeanValue = (i < m_FirstNumberOfComponents) ? m_FirstMeanValues[i] : m_WorkMeanValues[i-m_FirstNumberOfComponents];
    double referencePrecisionValue = (i < m_FirstNumberOfComponents) ? m_FirstPrecisionValues[i] : m_SecondPrecisionValues[i-m_FirstNumberOfComponents];
    double referenceMixingValue = ((i < m_FirstNumberOfComponents) ? m_FirstMixingValues[i] : m_SecondMixingValues[i-m_FirstNumberOfComponents]) / 2.0;

    double  sldIntegral = 0.0;
    double   ldIntegral = 0.0;
    double  jldIntegral = 0.0;
    double sldjIntegral = 0.0;
    double  ldjIntegral = 0.0;
    double    jIntegral = 0.0;

    for (unsigned int j = 0;j < numPoints;++j)
    {
      double quadWeight = m_QuadratureWeights[j];
      double quadPoint = referenceMeanValue + std::sqrt(2.0 / referencePrecisionValue) * m_QuadraturePoints[j];

      double logValue1 = m_FirstMixture.GetLogDensity(quadPoint);
      double logValue2 = m_WorkMixture.GetLogDensity(quadPoint);
      double logDifference = logValue1 - logValue2;
      double shiftLDJacobian = m_WorkMixture.GetLogDensityShiftDerivative(quadPoint);

      sldIntegral += quadWeight * logDifference * logDifference;
      ldIntegral += quadWeight * logDifference;
      jldIntegral += quadWeight * logDifference * shiftLDJacobian;
      if (i >= m_FirstNumberOfComponents)
      {
        double workScalar = quadWeight * logDifference * (quadPoint - referenceMeanValue);
        sldjIntegral += logDifference * workScalar;
        ldjIntegral += workScalar;
      }
      jIntegral += quadWeight * shiftLDJacobian;
    }

    totalSLDIntegral += referenceMixingValue * sldIntegral;
    totalLDIntegral += referenceMixingValue * ldIntegral;
    totalJLDIntegral += referenceMixingValue * jldIntegral;
    if (i >= m_FirstNumberOfComponents)
    {
      totalSLDJIntegral += referencePrecisionValue * referenceMixingValue * sldjIntegral;
      totalLDJIntegral += referencePrecisionValue * referenceMixingValue * ldjIntegral;
    }
    totalJIntegral += referenceMixingValue * jIntegral;
  }

  double costValue = totalSLDIntegral / std::sqrt(M_PI) - totalLDIntegral * totalLDIntegral / M_PI;

  m_ShiftDerivative = -2.0 * totalJLDIntegral;
  m_ShiftDerivative += totalSLDJIntegral;
  m_ShiftDerivative /= std::sqrt(M_PI);
  m_ShiftDerivative -= 2.0 * totalLDIntegral * (totalLDJIntegral - totalJIntegral) / M_PI;

  return costValue;
}
