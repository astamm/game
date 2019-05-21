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

  // Rcpp::Rcout << x << std::endl;
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

  m_LDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);
  m_SLDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);
  m_LDFDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);
  m_SLDFDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);
  m_LDSDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);
  m_SLDSDFunction.SetReferenceMixture(m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues, m_TotalNumberOfComponents);

  return numComponents;
}

// double SquaredBayesDistance::EvaluateSquaredDistance()
// {
//   m_LDFunction.SetInput1(m_FirstMixtureCopy);
//   m_LDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalLogIntegral = integrate(m_LDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   m_SLDFunction.SetInput1(m_FirstMixtureCopy);
//   m_SLDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalSquareLogIntegral = integrate(m_SLDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   m_LDFDFunction.SetInput1(m_FirstMixtureCopy);
//   m_LDFDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalLogIntegralWRTFirst = integrate(m_LDFDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   m_SLDFDFunction.SetInput1(m_FirstMixtureCopy);
//   m_SLDFDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalSquareLogIntegralWRTFirst = integrate(m_SLDFDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   m_LDSDFunction.SetInput1(m_FirstMixtureCopy);
//   m_LDSDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalLogIntegralWRTSecond = integrate(m_LDSDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   m_SLDSDFunction.SetInput1(m_FirstMixtureCopy);
//   m_SLDSDFunction.SetInput2(m_SecondMixtureCopy);
//   const double totalSquareLogIntegralWRTSecond = integrate(m_SLDSDFunction, m_LowerBound, m_UpperBound, m_ErrorEstimate, m_ErrorCode, m_NumberOfSubdivisions, m_AbsoluteTolerance, m_RelativeTolerance, m_QuadratureRule);
//
//   double costValue = totalSquareLogIntegral - totalLogIntegral * totalLogIntegral;
//
//   m_FirstShiftDerivative = totalSquareLogIntegralWRTFirst;
//   m_FirstShiftDerivative -= totalLogIntegralWRTFirst * totalLogIntegral;
//   m_FirstShiftDerivative *= 2.0;
//
//   m_SecondShiftDerivative = totalSquareLogIntegralWRTSecond;
//   m_SecondShiftDerivative -= totalLogIntegralWRTSecond * totalLogIntegral;
//   m_SecondShiftDerivative *= 2.0;
//
//   return costValue;
// }

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
      double weightValue = m_QuadratureWeights[j];
      // if (weightValue * referenceMixingValue * numPoints * m_TotalNumberOfComponents < 1.0e-6)
      //   continue;

      double inputValue = m_ReferenceMeanValues[i] + std::sqrt(2.0 / m_ReferencePrecisionValues[i]) * m_QuadraturePoints[j];

      // if (std::abs(inputValue) > 2500.0)
      //   continue;

      // double logValue1 = m_FirstMixtureCopy.GetLogDensity(inputValue);
      // double logValue2 = m_SecondMixtureCopy.GetLogDensity(inputValue);
      // double logDifference = logValue1 - logValue2;
      double logDifference = m_FirstMixtureCopy.GetLogDifference(inputValue, m_SecondMixtureCopy);

      // if (logValue1 < 1.0e-16 || logValue2 < 1.0e-16)
      //   continue;

      // Rcpp::Rcout << "After log: " << logValue1 << " " << logValue2 << std::endl;

      // if (!std::isfinite(logValue1) || !std::isfinite(logValue2))
      //   Rcpp::stop("NAN");

      logIntegral += weightValue * logDifference;
      squareLogIntegral += weightValue * logDifference * logDifference;

      double firstJacobianValue = m_FirstMixtureCopy.GetLogDensityShiftDerivative(inputValue);
      logIntegralWRTFirst += weightValue * firstJacobianValue;
      squareLogIntegralWRTFirst += weightValue * firstJacobianValue * logDifference;

      double secondJacobianValue = -m_SecondMixtureCopy.GetLogDensityShiftDerivative(inputValue);
      logIntegralWRTSecond += weightValue * secondJacobianValue;
      squareLogIntegralWRTSecond += weightValue * secondJacobianValue * logDifference;
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
