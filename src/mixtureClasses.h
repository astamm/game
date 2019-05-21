#pragma once

#include <Rcpp.h>

class GenericMixture
{
public:
  typedef std::vector<double> VectorType;

  virtual void SetInput(const Rcpp::DataFrame &inputData) = 0;
  virtual double GetLogDifference(const double inputValue, const GenericMixture &rhs) const = 0;
  virtual double GetLogDensity(const double inputValue) const = 0;
  virtual double GetLogDensityShiftDerivative(const double inputValue) const = 0;

  void SetNumberOfComponents(const unsigned int num) {m_NumberOfComponents = num;}
  double GetNumberOfComponents() const {return m_NumberOfComponents;}

  void SetMeanValues(const VectorType &x) {m_MeanValues = x;}
  VectorType GetMeanValues() const {return m_MeanValues;}

  void SetPrecisionValues(const VectorType &x) {m_PrecisionValues = x;}
  VectorType GetPrecisionValues() const {return m_PrecisionValues;}

  void SetMixingValues(const VectorType &x) {m_MixingValues = x;}
  VectorType GetMixingValues() const {return m_MixingValues;}

  double GetMean();

  GenericMixture() {}
  virtual ~GenericMixture() {}

protected:
  VectorType m_MeanValues, m_PrecisionValues, m_MixingValues;
  unsigned int m_NumberOfComponents;
};

class GaussianMixture: public GenericMixture
{
public:
  void SetInput(const Rcpp::DataFrame &inputData);
  virtual double GetLogDifference(const double inputValue, const GenericMixture &rhs) const;
  double GetLogDensity(const double inputValue) const;
  double GetLogDensityShiftDerivative(const double inputValue) const;

private:
  mutable VectorType m_WorkVector;
};
