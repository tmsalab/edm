#include "ecdm.h"

//' Compute the Odds Ratio
//'
//' Calculates the Odds Ratio
//' @param N  Number of Observations
//' @param J  Number of Assessment Items
//' @param Yt Estimated binary responses to assessements in \code{matrix} form
//'           with dimensions \eqn{N \times J}{N x J}.
//' @return A `matrix` with dimensions \eqn{J \times J}{J x J}.
//' @export
// [[Rcpp::export]]
arma::mat OddsRatio(unsigned int N, unsigned int J, const arma::mat &Yt)
{
    arma::mat M2_temp = arma::zeros<arma::mat>(J, J);
    for (unsigned int j1 = 0; j1 < J - 1; ++j1) {
        for (unsigned int j2 = j1 + 1; j2 < J; ++j2) {
            double n11 = arma::accu(Yt.col(j1) % Yt.col(j2));
            double n00 = arma::accu((1. - Yt.col(j1)) % (1. - Yt.col(j2)));
            double n10 = arma::accu(Yt.col(j1) % (1. - Yt.col(j2)));
            double n01 = N - n11 - n00 - n10;
            M2_temp(j1, j2) = (n11 * n00) / (n10 * n01);
        }
    }
    return M2_temp;
}

//' Bijection Vector
//'
//' @param K number of levels
//' @return A \code{vec} with length \eqn{K}.
//' @export
// [[Rcpp::export]]
arma::vec bijectionvector(unsigned int K)
{
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
arma::vec inv_bijectionvector(unsigned int K, double CL)
{
    arma::vec alpha(K);

    for (unsigned int k = 0; k < K; ++k) {
        double twopow = pow(2, K - k - 1);
        alpha(k) = (twopow <= CL);
        CL = CL - twopow * alpha(k);
    }

    return alpha;
}
