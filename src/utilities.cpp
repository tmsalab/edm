#include "ecdm.h"

//' Bijection Vector
//'
//' @param K number of levels
//' @return A \code{vec} with length \eqn{K}.
//' @export
// [[Rcpp::export]]
arma::vec bijectionvector(unsigned int K) {
    arma::vec vv(K);

    for (unsigned int k = 0; k < K; ++k) {
        vv(k) = pow(2, K - k - 1);
    }

    return vv;
}

//' Inverse Bijection Vector
//'
//' @param CL A \code{double} that controls ...
//' @inheritParams bijectionvector
//' @return A \code{vec} with length \eqn{K}.
//' @export
// [[Rcpp::export]]
arma::vec inv_bijectionvector(unsigned int K, double CL) {
    arma::vec alpha(K);

    for (unsigned int k = 0; k < K; ++k) {
        double twopow = pow(2, K - k - 1);
        alpha(k) = (twopow <= CL);
        CL = CL - twopow * alpha(k);
    }

    return alpha;
}
