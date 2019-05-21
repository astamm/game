#include "mixtureClasses.h"

double GenericMixture::GetMean()
{
  double meanValue = 0.0;

  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
    meanValue += m_MixingValues[i] * m_MeanValues[i];

  return meanValue;
}

void GaussianMixture::SetInput(const Rcpp::DataFrame &inputData)
{
  this->SetNumberOfComponents(inputData.nrows());
  this->SetMeanValues(Rcpp::as<VectorType>(inputData["mean"]));
  this->SetPrecisionValues(Rcpp::as<VectorType>(inputData["precision"]));
  this->SetMixingValues(Rcpp::as<VectorType>(inputData["mixing"]));
}

double GaussianMixture::GetLogDifference(const double inputValue, const GenericMixture &rhs) const
{
  m_WorkVector.resize(m_NumberOfComponents + rhs.GetNumberOfComponents());

  unsigned int indexOfMaximalExponent = 0;
  double maximalExponent = 0.0;

  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
  {
    double workValue = -m_PrecisionValues[i] * (inputValue - m_MeanValues[i]) * (inputValue - m_MeanValues[i]) / 2.0;

    // if (workValue < maximalExponent || i == 0)
    // {
    //   maximalExponent = workValue;
    //   indexOfMaximalExponent = i;
    // }

    m_WorkVector[i] = workValue;
  }

  for (unsigned int i = 0;i < rhs.GetNumberOfComponents();++i)
  {
    double workValue = -rhs.GetPrecisionValues()[i] * (inputValue - rhs.GetMeanValues()[i]) * (inputValue - rhs.GetMeanValues()[i]) / 2.0;

    if (workValue > maximalExponent || i == 0)
    {
      maximalExponent = workValue;
      indexOfMaximalExponent = i + m_NumberOfComponents;
    }

    m_WorkVector[i + m_NumberOfComponents] = workValue;
  }

  double numerValue = -1.0;
  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
   numerValue += m_MixingValues[i] * std::sqrt(m_PrecisionValues[i] / (2.0 * M_PI)) * std::exp(m_WorkVector[i] - maximalExponent);

  double denomValue = 1.0;
  for (unsigned int i = 0;i < rhs.GetNumberOfComponents();++i)
  {
    if (i + m_NumberOfComponents == indexOfMaximalExponent)
      continue;
    double workValue = rhs.GetMixingValues()[i] * std::sqrt(rhs.GetPrecisionValues()[i] / (2.0 * M_PI)) * std::exp(m_WorkVector[i + m_NumberOfComponents] - maximalExponent);
    numerValue -= workValue;
    denomValue += workValue;
  }


  return std::log1p(numerValue / denomValue);
}

double GaussianMixture::GetLogDensity(const double inputValue) const
{
  m_WorkVector.resize(m_NumberOfComponents);

  unsigned int indexOfMaximalExponent = 0;
  double maximalExponent = 0.0;
  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
  {
    double workValue = -m_PrecisionValues[i] * (inputValue - m_MeanValues[i]) * (inputValue - m_MeanValues[i]) / 2.0;

    if (workValue > maximalExponent || i == 0)
    {
      maximalExponent = workValue;
      indexOfMaximalExponent = i;
    }

    if (workValue < -20.0)
      continue;

    m_WorkVector[i] = 0.5 * std::log(m_PrecisionValues[i] / (2.0 * M_PI)) + workValue;
  }

  double insideLogValue = 1.0;

  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
  {
    if (i == indexOfMaximalExponent)
      continue;

    insideLogValue += std::exp(m_WorkVector[i] - m_WorkVector[indexOfMaximalExponent]) * m_MixingValues[i] / m_MixingValues[indexOfMaximalExponent];
  }

  return std::log(m_MixingValues[indexOfMaximalExponent]) + m_WorkVector[indexOfMaximalExponent] + std::log(insideLogValue);
}

double GaussianMixture::GetLogDensityShiftDerivative(const double inputValue) const
{
  m_WorkVector.resize(m_NumberOfComponents);

  double maximalExponent = 0.0;
  unsigned int indexOfMaximalExponent = 0;

  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
  {
    double workScalar = -m_PrecisionValues[i] * (inputValue - m_MeanValues[i]) * (inputValue - m_MeanValues[i]) / 2.0;

    if (workScalar > maximalExponent || i == 0)
    {
      maximalExponent = workScalar;
      indexOfMaximalExponent = i;
    }

    m_WorkVector[i] = workScalar;
  }

  double sumWeights = 0.0;
  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
  {
    double workScalar = m_MixingValues[i] * std::sqrt(m_PrecisionValues[i] / (2.0 * M_PI)) * std::exp(m_WorkVector[i] - maximalExponent);
    m_WorkVector[i] = workScalar;
    sumWeights += workScalar;
  }

  if (sumWeights < 1.0e-16)
  {
    sumWeights = (double)m_NumberOfComponents;
    for (unsigned int i = 0;i < m_NumberOfComponents;++i)
      m_WorkVector[i] = 1.0;
  }

  double outputValue = 0.0;

  for (unsigned int i = 0;i < m_NumberOfComponents;++i)
    outputValue += m_PrecisionValues[i] * (inputValue - m_MeanValues[i]) * m_WorkVector[i];

  outputValue /= sumWeights;

  return outputValue;
}
