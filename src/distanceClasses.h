// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

class SquaredBayesDistance: public MFuncGrad
{
public:
  typedef std::vector<double> VectorType;
  typedef Eigen::VectorXd ParametersType;

  double f_grad(Constvec& x, Refvec grad);

  // Set observed mixtures and reference mixture
  void SetDataPoints(const Rcpp::List &inputData);
  void SetInput1(const unsigned int index) {this->SetInputValues(0, index);}
  void SetInput2(const unsigned int index) {this->SetInputValues(1, index);}
  void SetQuadraturePoints(const VectorType &points) {m_QuadraturePoints = points;}
  void SetQuadratureWeights(const VectorType &weights) {m_QuadratureWeights = weights;}

protected:
  void SetInputValues(const unsigned int input, const unsigned int index);
  unsigned int SetReferenceModel(const Rcpp::List &inputData);
  double EvaluateSquaredDistance(
      const VectorType &firstMeanValues, const VectorType &firstPrecisionValues, const VectorType &firstMixingValues,
      const VectorType &secondMeanValues, const VectorType &secondPrecisionValues, const VectorType &secondMixingValues,
      const VectorType &referenceMeanValues, const VectorType &referencePrecisionValues, const VectorType &referenceMixingValues,
      double &jacobianValueWRTFirst, double &jacobianValueWRTReference
  );
  double EvaluateLogDensityValue(
      const double inputValue,
      const VectorType &meanValues,
      const VectorType &precisionValues,
      const VectorType &mixingValues
  );
  double EvaluateDensityJacobianWRTMean(
      const double inputValue,
      const VectorType &meanValues,
      const VectorType &precisionValues,
      const VectorType &mixingValues
  );

private:
  VectorType m_FirstMeanValues, m_FirstPrecisionValues, m_FirstMixingValues;
  VectorType m_SecondMeanValues, m_SecondPrecisionValues, m_SecondMixingValues;
  VectorType m_QuadraturePoints, m_QuadratureWeights;
  VectorType m_ReferenceMeanValues, m_ReferencePrecisionValues, m_ReferenceMixingValues;
  VectorType m_WorkVector;
  Rcpp::DataFrame m_WorkTibble;
  Rcpp::List m_DataPoints;
  unsigned int m_TotalNumberOfComponents;
};
