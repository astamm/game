#pragma once

#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::NumericVector GetMixtureDensity(
                const Rcpp::NumericVector &inputValues,
                const Rcpp::DataFrame &inputModel,
                const bool logScale);

// [[Rcpp::export]]
Rcpp::NumericVector GetMean(
                const Rcpp::NumericVector &inputValues,
                const Rcpp::List &mixtureModels,
                const double alpha = 0.2);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistancesToMean(
                const Rcpp::List &inputData,
                const Rcpp::LogicalVector &subsetValues,
                const Rcpp::NumericVector &nodeValues,
                const Rcpp::NumericVector &weightValues);

// [[Rcpp::export]]
Rcpp::NumericVector GetSquaredDistanceMatrix(
                const Rcpp::List &inputData,
                const std::vector<double> &nodeValues,
                const std::vector<double> &weightValues);

// [[Rcpp::export]]
double GetMeanRawMoment(
                const Rcpp::List &inputData,
                const unsigned int order,
                const Rcpp::NumericVector &nodeValues,
                const Rcpp::NumericVector &weightValues);

double GetLogDensityValue(
                const double inputValue,
                const Rcpp::NumericVector &meanValues,
                const Rcpp::NumericVector &precisionValues,
                const Rcpp::NumericVector &mixingValues,
                std::vector<double> &workVector);

// [[Rcpp::export]]
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
                const unsigned int numComponents,
                const Rcpp::NumericVector &nodeValues,
                const Rcpp::NumericVector &weightValues,
                std::vector<double> &workVector);

unsigned int BuildReferenceModel(
                const Rcpp::List &inputData,
                std::vector<double> &referenceMeanValues,
                std::vector<double> &referencePrecisionValues,
                std::vector<double> &referenceMixingValues);
