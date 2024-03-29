// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// NN_OCN
List NN_OCN(int dimX, int dimY, bool periodicBoundaries, IntegerMatrix movement);
RcppExport SEXP _OCNet_NN_OCN(SEXP dimXSEXP, SEXP dimYSEXP, SEXP periodicBoundariesSEXP, SEXP movementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dimX(dimXSEXP);
    Rcpp::traits::input_parameter< int >::type dimY(dimYSEXP);
    Rcpp::traits::input_parameter< bool >::type periodicBoundaries(periodicBoundariesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type movement(movementSEXP);
    rcpp_result_gen = Rcpp::wrap(NN_OCN(dimX, dimY, periodicBoundaries, movement));
    return rcpp_result_gen;
END_RCPP
}
// NN_river
List NN_river(int dimX, int dimY, bool periodicBoundaries, IntegerMatrix movement, IntegerVector toDEM, int nNodes);
RcppExport SEXP _OCNet_NN_river(SEXP dimXSEXP, SEXP dimYSEXP, SEXP periodicBoundariesSEXP, SEXP movementSEXP, SEXP toDEMSEXP, SEXP nNodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dimX(dimXSEXP);
    Rcpp::traits::input_parameter< int >::type dimY(dimYSEXP);
    Rcpp::traits::input_parameter< bool >::type periodicBoundaries(periodicBoundariesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type movement(movementSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type toDEM(toDEMSEXP);
    Rcpp::traits::input_parameter< int >::type nNodes(nNodesSEXP);
    rcpp_result_gen = Rcpp::wrap(NN_river(dimX, dimY, periodicBoundaries, movement, toDEM, nNodes));
    return rcpp_result_gen;
END_RCPP
}
// NN_FD
List NN_FD(int nNodes, int dimX, int dimY, List NeighbouringNodes, IntegerVector toDEM);
RcppExport SEXP _OCNet_NN_FD(SEXP nNodesSEXP, SEXP dimXSEXP, SEXP dimYSEXP, SEXP NeighbouringNodesSEXP, SEXP toDEMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nNodes(nNodesSEXP);
    Rcpp::traits::input_parameter< int >::type dimX(dimXSEXP);
    Rcpp::traits::input_parameter< int >::type dimY(dimYSEXP);
    Rcpp::traits::input_parameter< List >::type NeighbouringNodes(NeighbouringNodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type toDEM(toDEMSEXP);
    rcpp_result_gen = Rcpp::wrap(NN_FD(nNodes, dimX, dimY, NeighbouringNodes, toDEM));
    return rcpp_result_gen;
END_RCPP
}
// WSC
List WSC(int nNodes, List SC_to_FD, IntegerVector FD_to_SC, List NeighbouringNodes);
RcppExport SEXP _OCNet_WSC(SEXP nNodesSEXP, SEXP SC_to_FDSEXP, SEXP FD_to_SCSEXP, SEXP NeighbouringNodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nNodes(nNodesSEXP);
    Rcpp::traits::input_parameter< List >::type SC_to_FD(SC_to_FDSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type FD_to_SC(FD_to_SCSEXP);
    Rcpp::traits::input_parameter< List >::type NeighbouringNodes(NeighbouringNodesSEXP);
    rcpp_result_gen = Rcpp::wrap(WSC(nNodes, SC_to_FD, FD_to_SC, NeighbouringNodes));
    return rcpp_result_gen;
END_RCPP
}
// continue_FD_SC
List continue_FD_SC(IntegerVector IndexHeadpixel, IntegerVector FD_to_SC, List SC_to_FD, IntegerVector downNode);
RcppExport SEXP _OCNet_continue_FD_SC(SEXP IndexHeadpixelSEXP, SEXP FD_to_SCSEXP, SEXP SC_to_FDSEXP, SEXP downNodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type IndexHeadpixel(IndexHeadpixelSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type FD_to_SC(FD_to_SCSEXP);
    Rcpp::traits::input_parameter< List >::type SC_to_FD(SC_to_FDSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type downNode(downNodeSEXP);
    rcpp_result_gen = Rcpp::wrap(continue_FD_SC(IndexHeadpixel, FD_to_SC, SC_to_FD, downNode));
    return rcpp_result_gen;
END_RCPP
}
// paths_cpp
List paths_cpp(S4 OCN, IntegerVector whichNodes, String str, bool includePaths, bool includeDownstreamNode, bool includeUnconnectedPaths);
RcppExport SEXP _OCNet_paths_cpp(SEXP OCNSEXP, SEXP whichNodesSEXP, SEXP strSEXP, SEXP includePathsSEXP, SEXP includeDownstreamNodeSEXP, SEXP includeUnconnectedPathsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type OCN(OCNSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type whichNodes(whichNodesSEXP);
    Rcpp::traits::input_parameter< String >::type str(strSEXP);
    Rcpp::traits::input_parameter< bool >::type includePaths(includePathsSEXP);
    Rcpp::traits::input_parameter< bool >::type includeDownstreamNode(includeDownstreamNodeSEXP);
    Rcpp::traits::input_parameter< bool >::type includeUnconnectedPaths(includeUnconnectedPathsSEXP);
    rcpp_result_gen = Rcpp::wrap(paths_cpp(OCN, whichNodes, str, includePaths, includeDownstreamNode, includeUnconnectedPaths));
    return rcpp_result_gen;
END_RCPP
}
