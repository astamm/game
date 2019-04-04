#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector GetLogDensityWRTLebesgue(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues);

// [[Rcpp::export]]
Rcpp::NumericVector GetLogDensityWRTGaussian(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    const double referenceMean,
    const double referencePrecision);

// [[Rcpp::export]]
Rcpp::NumericVector GetCenteredLogRatio(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    const double referenceMean,
    const double referencePrecision,
    const double centeringValue);

// [[Rcpp::export]]
Rcpp::NumericVector GetMean(
    const Rcpp::NumericVector &inputValues,
    const Rcpp::List &mixtureModels,
    const double alpha = 0.2);

// [[Rcpp::export]]
Rcpp::NumericVector GetMeanSquaredNorms(
    const Rcpp::List &mixtureModels,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistancesToMean(
    const Rcpp::List &mixtureModels,
    const Rcpp::LogicalVector &subsetValues,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistanceMatrix(
    const Rcpp::List &inputModels,
    const Rcpp::DataFrame &referenceModel,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues);

double GetLogDensityValue(
    const double inputValue,
    const Rcpp::NumericVector &meanValues,
    const Rcpp::NumericVector &precisionValues,
    const Rcpp::NumericVector &mixingValues,
    std::vector<double> &logDensityValues);

// [[Rcpp::export]]
Rcpp::List PerformReassignment(
    const Rcpp::IntegerVector &inputMemberships,
    const Rcpp::List &inputModels,
    const Rcpp::IntegerVector &referenceModels,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues,
    const double alpha = 0.2);

// [[Rcpp::export]]
Rcpp::IntegerVector ComputeNewReferences(
    const Rcpp::IntegerVector &inputReferences,
    const Rcpp::List &inputModels,
    const Rcpp::IntegerVector &membershipVector,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues);

double GetDistanceToMean(
    const Rcpp::List &mixtureModels,
    const unsigned int observationIndex,
    const Rcpp::IntegerVector &membershipValues,
    const unsigned int clusterIndex,
    const unsigned int referenceIndex,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues,
    const double alpha = 0.2);

// [[Rcpp::export]]
double GetDistance(
    const Rcpp::DataFrame &firstModel,
    const Rcpp::DataFrame &secondModel,
    const Rcpp::DataFrame &referenceModel,
    const Rcpp::NumericVector &nodeValues,
    const Rcpp::NumericVector &weightValues);
