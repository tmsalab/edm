#ifndef ECDM_UTILITIES_H
#define ECDM_UTILITIES_H
arma::vec bijectionvector(unsigned int K);
arma::vec inv_bijectionvector(unsigned int K, double CL);
arma::mat OddsRatio(unsigned int N, unsigned int J, const arma::mat &Yt);
#endif
