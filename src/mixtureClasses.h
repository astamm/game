#include <Rcpp.h>

class GenericMixture
{
public:
  typedef std::vector<double> VectorType;

  virtual void SetInput(const Rcpp::DataFrame &inputData) = 0;
  virtual double GetLogDensity(const double inputValue) = 0;
  virtual double GetLogDensityShiftDerivative(const double inputValue) = 0;

  void SetNumberOfComponents(const unsigned int num) {m_NumberOfComponents = num;}
  double GetNumberOfComponents() {return m_NumberOfComponents;}

  void SetMeanValues(const VectorType &x) {m_MeanValues = x;}
  VectorType GetMeanValues() {return m_MeanValues;}

  void SetPrecisionValues(const VectorType &x) {m_PrecisionValues = x;}
  VectorType GetPrecisionValues() {return m_PrecisionValues;}

  void SetMixingValues(const VectorType &x) {m_MixingValues = x;}
  VectorType GetMixingValues() {return m_MixingValues;}

  double GetMean();

protected:
  VectorType m_MeanValues, m_PrecisionValues, m_MixingValues;
  unsigned int m_NumberOfComponents;
};

class GaussianMixture: public GenericMixture
{
public:
  void SetInput(const Rcpp::DataFrame &inputData);
  double GetLogDensity(const double inputValue);
  double GetLogDensityShiftDerivative(const double inputValue);

private:
  VectorType m_WorkVector;
};
