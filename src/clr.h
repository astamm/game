#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector GetMean(
        const Rcpp::NumericVector &inputValues,
        const Rcpp::List &mixtureModels,
        const double alpha = 0.2);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistancesToMean(
        const Rcpp::List &mixtureModels,
        const Rcpp::LogicalVector &subsetValues,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistanceMatrix(
        const Rcpp::List &inputData,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
double GetMeanNormalizationFactor(
        const Rcpp::List &inputData,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
double GetMeanExpectedValue(
        const Rcpp::List &inputData,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues);

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

double GetLogDensityValue(
        const double inputValue,
        const Rcpp::NumericVector &meanValues,
        const Rcpp::NumericVector &precisionValues,
        const Rcpp::NumericVector &mixingValues,
        std::vector<double> &workVector);

double GetDistanceToMean(
        const Rcpp::List &mixtureModels,
        const unsigned int observationIndex,
        const Rcpp::IntegerVector &membershipValues,
        const unsigned int clusterIndex,
        const unsigned int referenceIndex,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues,
        const double alpha = 0.2);

double GetSquaredDistance(
        const Rcpp::NumericVector &firstMeanValues,
        const Rcpp::NumericVector &firstPrecisionValues,
        const Rcpp::NumericVector &firstMixingValues,
        const Rcpp::NumericVector &secondMeanValues,
        const Rcpp::NumericVector &secondPrecisionValues,
        const Rcpp::NumericVector &secondMixingValues,
        const std::vector<double> &referenceMeanValues,
        const std::vector<double> &referencePrecisionValues,
        const std::vector<double> &referenceMixingValues,
        const Rcpp::NumericVector &nodeValues,
        const Rcpp::NumericVector &weightValues,
        std::vector<double> &workVector);

unsigned int BuildReferenceModel(
        const Rcpp::List &inputData,
        std::vector<double> &referenceMeanValues,
        std::vector<double> &referencePrecisionValues,
        std::vector<double> &referenceMixingValues);
