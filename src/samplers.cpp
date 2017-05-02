#include "ecdm.h"

//' @title Generate Multinomial Random Variable
//' @description Sample a multinomial random variable for given probabilities.
//' @usage rmultinomial(ps)
//' @param ps A \code{vector} for the probability of each category.
//' @return A \code{vector} from a multinomial with probability ps.
//' @author Steven Andrew Culpepper
// [[Rcpp::export]]
double rmultinomial(const arma::vec &ps) {
    unsigned int C = ps.n_elem;
    double u = R::runif(0, 1);
    arma::vec cps = cumsum(ps);
    arma::vec Ips = arma::zeros<arma::vec>(C);
    Ips.elem(arma::find(cps < u)).fill(1.0);

    return sum(Ips);
}

//' @title Generate Dirichlet Random Variable
//' @description Sample a Dirichlet random variable.
//' @usage rDirichlet(deltas)
//' @param deltas A \code{vector} of Dirichlet parameters.
//' @return A \code{vector} from a Dirichlet.
//' @author Steven Andrew Culpepper
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec &deltas) {
    unsigned int C = deltas.n_elem;
    arma::vec Xgamma(C);

    // generating gamma(deltac,1)
    for (unsigned int c = 0; c < C; c++) {
        Xgamma(c) = R::rgamma(deltas(c), 1.0);
    }
    return Xgamma / sum(Xgamma);
}
