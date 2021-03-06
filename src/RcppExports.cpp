// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fa2_distance
double fa2_distance(NumericVector& n1, NumericVector& n2, double norm = 2.0);
RcppExport SEXP RforceAtlas2_fa2_distance(SEXP n1SEXP, SEXP n2SEXP, SEXP normSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector& >::type n1(n1SEXP );
        Rcpp::traits::input_parameter< NumericVector& >::type n2(n2SEXP );
        Rcpp::traits::input_parameter< double >::type norm(normSEXP );
        double __result = fa2_distance(n1, n2, norm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_getids
IntegerVector fa2_getids(IntegerMatrix& G);
RcppExport SEXP RforceAtlas2_fa2_getids(SEXP GSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix& >::type G(GSEXP );
        IntegerVector __result = fa2_getids(G);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_degree
IntegerMatrix fa2_degree(IntegerMatrix& G);
RcppExport SEXP RforceAtlas2_fa2_degree(SEXP GSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix& >::type G(GSEXP );
        IntegerMatrix __result = fa2_degree(G);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_degree_index
IntegerVector fa2_degree_index(IntegerMatrix& G);
RcppExport SEXP RforceAtlas2_fa2_degree_index(SEXP GSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix& >::type G(GSEXP );
        IntegerVector __result = fa2_degree_index(G);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_get_center
NumericVector fa2_get_center(NumericMatrix& pos);
RcppExport SEXP RforceAtlas2_fa2_get_center(SEXP posSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix& >::type pos(posSEXP );
        NumericVector __result = fa2_get_center(pos);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_forcei
double fa2_forcei(NumericVector n1, NumericVector n2, int d1, int d2, double s1 = 1.0, double s2 = 1.0, double kr = 1.0, double kr2 = 100.0, double ka = 1.0, bool nooverlap = false);
RcppExport SEXP RforceAtlas2_fa2_forcei(SEXP n1SEXP, SEXP n2SEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP krSEXP, SEXP kr2SEXP, SEXP kaSEXP, SEXP nooverlapSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type n1(n1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type n2(n2SEXP );
        Rcpp::traits::input_parameter< int >::type d1(d1SEXP );
        Rcpp::traits::input_parameter< int >::type d2(d2SEXP );
        Rcpp::traits::input_parameter< double >::type s1(s1SEXP );
        Rcpp::traits::input_parameter< double >::type s2(s2SEXP );
        Rcpp::traits::input_parameter< double >::type kr(krSEXP );
        Rcpp::traits::input_parameter< double >::type kr2(kr2SEXP );
        Rcpp::traits::input_parameter< double >::type ka(kaSEXP );
        Rcpp::traits::input_parameter< bool >::type nooverlap(nooverlapSEXP );
        double __result = fa2_forcei(n1, n2, d1, d2, s1, s2, kr, kr2, ka, nooverlap);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fa2_force
NumericMatrix fa2_force(NumericMatrix& pos, IntegerMatrix& G, NumericVector& size, double kr = 1.0, double kr2 = 100.0, double ka = 1.0, bool nooverlap = false);
RcppExport SEXP RforceAtlas2_fa2_force(SEXP posSEXP, SEXP GSEXP, SEXP sizeSEXP, SEXP krSEXP, SEXP kr2SEXP, SEXP kaSEXP, SEXP nooverlapSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix& >::type pos(posSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix& >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector& >::type size(sizeSEXP );
        Rcpp::traits::input_parameter< double >::type kr(krSEXP );
        Rcpp::traits::input_parameter< double >::type kr2(kr2SEXP );
        Rcpp::traits::input_parameter< double >::type ka(kaSEXP );
        Rcpp::traits::input_parameter< bool >::type nooverlap(nooverlapSEXP );
        NumericMatrix __result = fa2_force(pos, G, size, kr, kr2, ka, nooverlap);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
