// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include "mixtureClasses.h"

class GenericIntegrand: public Numer::Func
{
public:
  typedef GenericMixture::VectorType VectorType;

  void SetInput1(const GaussianMixture &input) {m_FirstMixture = input;}
  void SetInput2(const GaussianMixture &input) {m_SecondMixture = input;}
  void SetReferenceMixture(
      const VectorType &means,
      const VectorType &precisions,
      const VectorType &mixings,
      const unsigned int ncomponents
  );
  double operator() (const double& x) const = 0;

protected:
  GaussianMixture m_FirstMixture;
  GaussianMixture m_SecondMixture;
  GaussianMixture m_ReferenceMixture;
};

class LogDifference: public GenericIntegrand
{
public:
  double operator() (const double& x) const;
};

class SquaredLogDifference: public GenericIntegrand
{
public:
  double operator() (const double& x) const;
};

class LogDifferenceFirstDerivative: public GenericIntegrand
{
public:
  double operator() (const double& x) const;
};

class SquaredLogDifferenceFirstDerivative: public GenericIntegrand
{
public:
  double operator() (const double& x) const;
};

class LogDifferenceSecondDerivative: public GenericIntegrand
{
public:
  double operator()(const double& x) const;
};

class SquaredLogDifferenceSecondDerivative: public GenericIntegrand
{
public:
  double operator() (const double& x) const;
};
