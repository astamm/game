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

  void SetInput1(const Rcpp::DataFrame &x);
  void SetInput2(const Rcpp::DataFrame &x);
  void SetQuadraturePoints(const VectorType &points) {m_QuadraturePoints = points;}
  void SetQuadratureWeights(const VectorType &weights) {m_QuadratureWeights = weights;}

  double GetFirstMeanValue() {return m_FirstMixture.GetMean();}
  double GetSecondMeanValue() {return m_SecondMixture.GetMean();}

protected:
  double EvaluateSquaredDistance();

private:
  GaussianMixture m_FirstMixture, m_SecondMixture;
  GaussianMixture m_SecondMixtureCopy;
  VectorType m_QuadraturePoints, m_QuadratureWeights;
  VectorType m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues;
  double m_SecondShiftDerivative;
  VectorType m_WorkVector;
  Rcpp::DataFrame m_WorkTibble;
  GaussianMixture m_WorkMixture;
  unsigned int m_ReferenceNumberOfComponents;
};
