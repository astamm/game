// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// GetLogDensityWRTLebesgue
Rcpp::NumericVector GetLogDensityWRTLebesgue(const Rcpp::NumericVector& inputValues, const Rcpp::NumericVector& meanValues, const Rcpp::NumericVector& precisionValues, const Rcpp::NumericVector& mixingValues);
RcppExport SEXP _game_GetLogDensityWRTLebesgue(SEXP inputValuesSEXP, SEXP meanValuesSEXP, SEXP precisionValuesSEXP, SEXP mixingValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputValues(inputValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type meanValues(meanValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type precisionValues(precisionValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mixingValues(mixingValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogDensityWRTLebesgue(inputValues, meanValues, precisionValues, mixingValues));
    return rcpp_result_gen;
END_RCPP
}
// GetLogDensityWRTGaussian
Rcpp::NumericVector GetLogDensityWRTGaussian(const Rcpp::NumericVector& inputValues, const Rcpp::NumericVector& meanValues, const Rcpp::NumericVector& precisionValues, const Rcpp::NumericVector& mixingValues, const double referenceMean, const double referencePrecision);
RcppExport SEXP _game_GetLogDensityWRTGaussian(SEXP inputValuesSEXP, SEXP meanValuesSEXP, SEXP precisionValuesSEXP, SEXP mixingValuesSEXP, SEXP referenceMeanSEXP, SEXP referencePrecisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputValues(inputValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type meanValues(meanValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type precisionValues(precisionValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mixingValues(mixingValuesSEXP);
    Rcpp::traits::input_parameter< const double >::type referenceMean(referenceMeanSEXP);
    Rcpp::traits::input_parameter< const double >::type referencePrecision(referencePrecisionSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogDensityWRTGaussian(inputValues, meanValues, precisionValues, mixingValues, referenceMean, referencePrecision));
    return rcpp_result_gen;
END_RCPP
}
// GetCenteredLogRatio
Rcpp::NumericVector GetCenteredLogRatio(const Rcpp::NumericVector& inputValues, const Rcpp::NumericVector& meanValues, const Rcpp::NumericVector& precisionValues, const Rcpp::NumericVector& mixingValues, const double referenceMean, const double referencePrecision, const double centeringValue);
RcppExport SEXP _game_GetCenteredLogRatio(SEXP inputValuesSEXP, SEXP meanValuesSEXP, SEXP precisionValuesSEXP, SEXP mixingValuesSEXP, SEXP referenceMeanSEXP, SEXP referencePrecisionSEXP, SEXP centeringValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputValues(inputValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type meanValues(meanValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type precisionValues(precisionValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mixingValues(mixingValuesSEXP);
    Rcpp::traits::input_parameter< const double >::type referenceMean(referenceMeanSEXP);
    Rcpp::traits::input_parameter< const double >::type referencePrecision(referencePrecisionSEXP);
    Rcpp::traits::input_parameter< const double >::type centeringValue(centeringValueSEXP);
    rcpp_result_gen = Rcpp::wrap(GetCenteredLogRatio(inputValues, meanValues, precisionValues, mixingValues, referenceMean, referencePrecision, centeringValue));
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
// GetMeanSquaredNorms
Rcpp::NumericVector GetMeanSquaredNorms(const Rcpp::List& mixtureModels, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetMeanSquaredNorms(SEXP mixtureModelsSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type mixtureModels(mixtureModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetMeanSquaredNorms(mixtureModels, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetSquaredDistancesToMean
Rcpp::NumericVector GetSquaredDistancesToMean(const Rcpp::List& mixtureModels, const Rcpp::LogicalVector& subsetValues, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetSquaredDistancesToMean(SEXP mixtureModelsSEXP, SEXP subsetValuesSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type mixtureModels(mixtureModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type subsetValues(subsetValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSquaredDistancesToMean(mixtureModels, subsetValues, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetSquaredDistanceMatrix
Rcpp::NumericVector GetSquaredDistanceMatrix(const Rcpp::List& inputModels, const Rcpp::DataFrame& referenceModel, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetSquaredDistanceMatrix(SEXP inputModelsSEXP, SEXP referenceModelSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputModels(inputModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type referenceModel(referenceModelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSquaredDistanceMatrix(inputModels, referenceModel, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// PerformReassignment
Rcpp::List PerformReassignment(const Rcpp::IntegerVector& inputMemberships, const Rcpp::List& inputModels, const Rcpp::IntegerVector& referenceModels, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues, const double alpha);
RcppExport SEXP _game_PerformReassignment(SEXP inputMembershipsSEXP, SEXP inputModelsSEXP, SEXP referenceModelsSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type inputMemberships(inputMembershipsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputModels(inputModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type referenceModels(referenceModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(PerformReassignment(inputMemberships, inputModels, referenceModels, nodeValues, weightValues, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ComputeNewReferences
Rcpp::IntegerVector ComputeNewReferences(const Rcpp::IntegerVector& inputReferences, const Rcpp::List& inputModels, const Rcpp::IntegerVector& membershipVector, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_ComputeNewReferences(SEXP inputReferencesSEXP, SEXP inputModelsSEXP, SEXP membershipVectorSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type inputReferences(inputReferencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type inputModels(inputModelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type membershipVector(membershipVectorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeNewReferences(inputReferences, inputModels, membershipVector, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}
// GetDistance
double GetDistance(const Rcpp::DataFrame& firstModel, const Rcpp::DataFrame& secondModel, const Rcpp::DataFrame& referenceModel, const Rcpp::NumericVector& nodeValues, const Rcpp::NumericVector& weightValues);
RcppExport SEXP _game_GetDistance(SEXP firstModelSEXP, SEXP secondModelSEXP, SEXP referenceModelSEXP, SEXP nodeValuesSEXP, SEXP weightValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type firstModel(firstModelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type secondModel(secondModelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type referenceModel(referenceModelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nodeValues(nodeValuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weightValues(weightValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetDistance(firstModel, secondModel, referenceModel, nodeValues, weightValues));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_game_GetLogDensityWRTLebesgue", (DL_FUNC) &_game_GetLogDensityWRTLebesgue, 4},
    {"_game_GetLogDensityWRTGaussian", (DL_FUNC) &_game_GetLogDensityWRTGaussian, 6},
    {"_game_GetCenteredLogRatio", (DL_FUNC) &_game_GetCenteredLogRatio, 7},
    {"_game_GetMean", (DL_FUNC) &_game_GetMean, 3},
    {"_game_GetMeanSquaredNorms", (DL_FUNC) &_game_GetMeanSquaredNorms, 3},
    {"_game_GetSquaredDistancesToMean", (DL_FUNC) &_game_GetSquaredDistancesToMean, 4},
    {"_game_GetSquaredDistanceMatrix", (DL_FUNC) &_game_GetSquaredDistanceMatrix, 4},
    {"_game_PerformReassignment", (DL_FUNC) &_game_PerformReassignment, 6},
    {"_game_ComputeNewReferences", (DL_FUNC) &_game_ComputeNewReferences, 5},
    {"_game_GetDistance", (DL_FUNC) &_game_GetDistance, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_game(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}