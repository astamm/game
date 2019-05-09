// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// GetMixtureDensity
Rcpp::NumericVector GetMixtureDensity(const Rcpp::NumericVector& inputValues, const Rcpp::DataFrame& inputModel, const bool logScale);
RcppExport SEXP _game_GetMixtureDensity(SEXP inputValuesSEXP, SEXP inputModelSEXP, SEXP logScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputValues(inputValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type inputModel(inputModelSEXP);
    Rcpp::traits::input_parameter< const bool >::type logScale(logScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(GetMixtureDensity(inputValues, inputModel, logScale));
    return rcpp_result_gen;
END_RCPP
}
// GetMean
Rcpp::NumericVector GetMean(const Rcpp::NumericVector& inputValues, const Rcpp::List& mixtureModels, const double alpha);
RcppExport SEXP _game_GetMean(SEXP inputValuesSEXP, SEXP mixtureModelsSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputValues(inputValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type mixtureModels(mixtureModelsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(GetMean(inputValues, mixtureModels, alpha));
    return rcpp_result_gen;
END_RCPP
}
// GetSquaredDistancesToMean
Rcpp::NumericVector GetSquaredDistancesToMean(const Rcpp::List& inputData, const Rcpp::LogicalVector& subsetValues, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetSquaredDistancesToMean(SEXP inputDataSEXP, SEXP subsetValuesSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputData(inputDataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type subsetValues(subsetValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSquaredDistancesToMean(inputData, subsetValues, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetSquaredDistanceMatrix
Rcpp::NumericVector GetSquaredDistanceMatrix(const Rcpp::List& inputData, const std::vector<double>& nodeValues, const std::vector<double>& weightValues);
RcppExport SEXP _game_GetSquaredDistanceMatrix(SEXP inputDataSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputData(inputDataSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSquaredDistanceMatrix(inputData, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetMeanRawMoment
double GetMeanRawMoment(const Rcpp::List& inputData, const unsigned int order, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetMeanRawMoment(SEXP inputDataSEXP, SEXP orderSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputData(inputDataSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetMeanRawMoment(inputData, order, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetSquaredDistance
double GetSquaredDistance(const Rcpp::NumericVector& firstMeanValues, const Rcpp::NumericVector& firstPrecisionValues, const Rcpp::NumericVector& firstMixingValues, const Rcpp::NumericVector& secondMeanValues, const Rcpp::NumericVector& secondPrecisionValues, const Rcpp::NumericVector& secondMixingValues, const std::vector<double>& referenceMeanValues, const std::vector<double>& referencePrecisionValues, const std::vector<double>& referenceMixingValues, const unsigned int numComponents, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues, std::vector<double>& workVector);
RcppExport SEXP _game_GetSquaredDistance(SEXP firstMeanValuesSEXP, SEXP firstPrecisionValuesSEXP, SEXP firstMixingValuesSEXP, SEXP secondMeanValuesSEXP, SEXP secondPrecisionValuesSEXP, SEXP secondMixingValuesSEXP, SEXP referenceMeanValuesSEXP, SEXP referencePrecisionValuesSEXP, SEXP referenceMixingValuesSEXP, SEXP numComponentsSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP, SEXP workVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type firstMeanValues(firstMeanValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type firstPrecisionValues(firstPrecisionValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type firstMixingValues(firstMixingValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type secondMeanValues(secondMeanValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type secondPrecisionValues(secondPrecisionValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type secondMixingValues(secondMixingValuesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type referenceMeanValues(referenceMeanValuesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type referencePrecisionValues(referencePrecisionValuesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type referenceMixingValues(referenceMixingValuesSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type numComponents(numComponentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type workVector(workVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSquaredDistance(firstMeanValues, firstPrecisionValues, firstMixingValues, secondMeanValues, secondPrecisionValues, secondMixingValues, referenceMeanValues, referencePrecisionValues, referenceMixingValues, numComponents, nodeValues, weightValues, workVector));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_game_GetMixtureDensity", (DL_FUNC) &_game_GetMixtureDensity, 3},
    {"_game_GetMean", (DL_FUNC) &_game_GetMean, 3},
    {"_game_GetSquaredDistancesToMean", (DL_FUNC) &_game_GetSquaredDistancesToMean, 4},
    {"_game_GetSquaredDistanceMatrix", (DL_FUNC) &_game_GetSquaredDistanceMatrix, 3},
    {"_game_GetMeanRawMoment", (DL_FUNC) &_game_GetMeanRawMoment, 4},
    {"_game_GetSquaredDistance", (DL_FUNC) &_game_GetSquaredDistance, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_game(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
