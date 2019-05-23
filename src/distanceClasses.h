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

  double GetMean1() {return m_FirstMixture.GetMean();}
  double GetMean2() {return m_SecondMixture.GetMean();}

protected:
  double EvaluateSquaredDistance();

private:
  GaussianMixture m_FirstMixture, m_SecondMixture;
  GaussianMixture m_WorkMixture;
  VectorType m_QuadraturePoints, m_QuadratureWeights;
  double m_SecondShiftDerivative;
  unsigned int m_FirstNumberOfComponents, m_SecondNumberOfComponents;
  VectorType m_FirstMeanValues, m_FirstPrecisionValues, m_FirstMixingValues;
  VectorType m_SecondMeanValues, m_SecondPrecisionValues, m_SecondMixingValues;
  VectorType m_WorkMeanValues;
};
