// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
#include "mixtureClasses.h"
#include "integrandClasses.h"

class SquaredBayesDistance: public Numer::MFuncGrad
{
public:
  typedef std::vector<double> VectorType;
  typedef Eigen::VectorXd ParametersType;

  double f_grad(Numer::Constvec& x, Numer::Refvec grad);

  // Set observed mixtures and reference mixture
  void SetDataPoints(const Rcpp::List &inputData);
  void SetInput1(const unsigned int index) {this->SetInputValues(0, index);}
  void SetInput2(const unsigned int index) {this->SetInputValues(1, index);}
  void SetQuadraturePoints(const VectorType &points) {m_QuadraturePoints = points;}
  void SetQuadratureWeights(const VectorType &weights) {m_QuadratureWeights = weights;}

  double GetFirstMeanValue() {return m_FirstMixture.GetMean();}
  double GetSecondMeanValue() {return m_SecondMixture.GetMean();}
  double GetReferenceMeanValue() {return m_ReferenceMeanValue;}

  void SetLowerBound(const double lb) {m_LowerBound = lb;}
  void SetUpperBound(const double ub) {m_UpperBound = ub;}

  SquaredBayesDistance() : Numer::MFuncGrad()
  {
    m_LowerBound = -2500.0;
    m_UpperBound =  2500.0;
  }

protected:
  void SetInputValues(const unsigned int input, const unsigned int index);
  unsigned int SetReferenceModel(const Rcpp::List &inputData);
  double EvaluateSquaredDistance();

private:
  GaussianMixture m_FirstMixture, m_SecondMixture;
  GaussianMixture m_FirstMixtureCopy, m_SecondMixtureCopy;
  VectorType m_QuadraturePoints, m_QuadratureWeights;
  VectorType m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues;
  double m_ReferenceMeanValue;
  double m_FirstShiftDerivative, m_SecondShiftDerivative;
  VectorType m_WorkVector;
  Rcpp::DataFrame m_WorkTibble;
  GaussianMixture m_WorkMixture;
  Rcpp::List m_DataPoints;
  unsigned int m_TotalNumberOfComponents;

  LogDifference m_LDFunction;
  SquaredLogDifference m_SLDFunction;
  LogDifferenceFirstDerivative m_LDFDFunction;
  SquaredLogDifferenceFirstDerivative m_SLDFDFunction;
  LogDifferenceSecondDerivative m_LDSDFunction;
  SquaredLogDifferenceSecondDerivative m_SLDSDFunction;

  double m_LowerBound, m_UpperBound;
  double m_ErrorEstimate;
  int m_ErrorCode;
  const int m_NumberOfSubdivisions = 100;
  const double m_AbsoluteTolerance = 1e-8;
  const double m_RelativeTolerance = 1e-6;
  const Numer::Integrator<double>::QuadratureRule m_QuadratureRule = Numer::Integrator<double>::GaussKronrod21;
};
