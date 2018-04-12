#include "ecdm.h"
#include <rgen.h>

// Remove after dealing with CLASS in place of alpha
arma::mat simrRUM_errum(unsigned int N, unsigned int J, unsigned int K,
                        const arma::mat &Q, const arma::mat &rstar,
                        const arma::vec &pistar, const arma::vec &CLASS)
{

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double kj;
    double aik;
    arma::vec pmu(J);
    arma::mat Y = arma::zeros<arma::mat>(N, J);
    for (unsigned int i = 0; i < N; i++) {
        arma::vec Yi = arma::zeros<arma::vec>(J);
        arma::vec pi = arma::ones<arma::vec>(J);
        arma::vec ui = arma::randu<arma::vec>(J);
        arma::vec alpha_i = inv_bijectionvector(K, CLASS(i));
        for (unsigned int j = 0; j < J; j++) {
            arma::uvec task_ij = find(Q.row(j) == 1);

            for (unsigned int k = 0; k < task_ij.n_elem; k++) {
                kj = task_ij(k);
                aik = alpha_i(kj);
                pi(j) = ((rstar(j, kj) * (1.0 - aik) + 1.0 * aik) * Q(j, kj) +
                         1.0 * (1.0 - Q(j, kj))) *
                        pi(j);
            }
            pi(j) = pistar(j) * pi(j);
        }
        pmu = pi - ui;
        Yi(arma::find(pmu > 0)).fill(1);
        Y.row(i) = Yi.t();
    }
    return Y;
}

arma::mat P_X_CL(unsigned int K, unsigned int nClass, const arma::vec &qj,
                 const arma::vec &gj, const arma::vec &sj)
{
    arma::mat Pmat = arma::ones<arma::mat>(nClass, nClass);
    arma::uvec qj_ones = find(qj == 1);
    arma::vec gj_ones = gj(qj_ones);
    arma::vec sj_ones = sj(qj_ones);
    double Kj = arma::accu(qj);
    for (unsigned int r = 0; r < nClass; r++) {
        arma::vec X_c = inv_bijectionvector(K, r);
        arma::vec X_c_ones = X_c(qj_ones);
        for (unsigned int c = 0; c < nClass; c++) {
            arma::vec alpha_c = inv_bijectionvector(K, c);
            arma::vec alpha_c_ones = alpha_c(qj_ones);
            for (unsigned int k = 0; k < Kj; k++) {
                double xk = X_c_ones(k);
                double ak = alpha_c_ones(k);
                double gk = gj_ones(k);
                double sk = sj_ones(k);
                Pmat(r, c) *= ((1. - sk) * ak + gk * (1. - ak)) * xk +
                              (sk * ak + (1. - gk) * (1. - ak)) * (1. - xk);
            }
        }
    }
    return Pmat;
}

// a single factor higher order version of of P_X_CL to also compute
// p(alpha|tau,lambda,xi)
double pis_HO(unsigned int K, unsigned int class_i, const arma::vec &tau,
              const arma::vec &lambda, const double xi_i,
              const arma::mat &alpha_class_tab)
{
    arma::vec alpha_c = alpha_class_tab.col(class_i);
    double palpha_c = 1.;
    for (unsigned int k = 0; k < K; k++) {
        double ak = alpha_c(k);
        palpha_c *= R::pnorm(0.0, tau(k) + lambda(k) * xi_i, 1., 1 - ak, 0);
    }
    return palpha_c;
}

arma::mat P_X_CL_Yeq0(unsigned int K, unsigned int nClass, const arma::vec &qj,
                      const arma::vec &gj, const arma::vec &sj)
{
    arma::mat Pmat = arma::ones<arma::mat>(nClass, nClass);
    arma::uvec qj_ones = find(qj == 1);
    double Kj = arma::accu(qj);
    for (unsigned int r = 0; r < nClass; r++) {
        arma::vec X_c = inv_bijectionvector(K, r);
        arma::vec X_c_ones = X_c(qj_ones);
        for (unsigned int c = 0; c < nClass; c++) {
            if (arma::accu(X_c % qj) == Kj) {
                Pmat(r, c) = .0;
            } else {
                arma::vec alpha_c = inv_bijectionvector(K, c);
                for (unsigned int k = 0; k < K; k++) {
                    double xk = X_c(k);
                    double ak = alpha_c(k);
                    double gk = gj(k);
                    double sk = sj(k);
                    double qk = qj(k);
                    Pmat(r, c) *=
                        (((1. - sk) * ak + gk * (1. - ak)) * xk +
                         (sk * ak + (1. - gk) * (1. - ak)) * (1. - xk)) *
                            qk +
                        (1. - qk) * .5;
                }
            }
        }
    }
    return Pmat;
}

arma::mat P_X_CL_Yeq1(unsigned int K, unsigned int nClass)
{
    arma::mat Pmat = arma::zeros<arma::mat>(nClass, nClass);
    for (unsigned int r = 0; r < nClass; r++) {
        arma::vec X_c = inv_bijectionvector(K, r);
        for (unsigned int c = 0; c < nClass; c++) {
            arma::vec qj = inv_bijectionvector(K, c);
            double Kj = arma::accu(qj);
            if (arma::accu(X_c % qj) == Kj) {
                Pmat(r, c) = 1. / pow(2, K - Kj);
            }
        }
    }
    return Pmat;
}

arma::cube X_CL_match(unsigned int K, unsigned int nClass)
{
    arma::cube match(nClass, nClass, K); // X by CL by match
    for (unsigned int r = 0; r < nClass; r++) {
        arma::vec X_c = inv_bijectionvector(K, r);
        for (unsigned int c = 0; c < nClass; c++) {
            arma::vec alpha_c = inv_bijectionvector(K, c);
            for (unsigned int k = 0; k < K; k++) {
                double xk = X_c(k);
                double ak = alpha_c(k);
                match(r, c, k) = 2 * xk + ak;
            }
        }
    }
    return match;
}

double rRUM_LLj_2(unsigned int N, unsigned int K, unsigned int J,
                  unsigned int nClass, const arma::vec &qj, const arma::vec &Yj,
                  const arma::vec &rstarj, double pistar, const arma::vec &pis)
{
    double LLj = 0;
    for (unsigned int i = 0; i < N; i++) {
        arma::vec pic(nClass);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            arma::vec ac = inv_bijectionvector(K, cc);
            double prod_ijk = 1.0;
            for (unsigned int k = 0; k < K; k++) {
                prod_ijk = pow(rstarj(k), qj(k) * (1. - ac(k))) * prod_ijk;
            }
            pic(cc) = pistar * prod_ijk * Yj(i) +
                      (1. - pistar * prod_ijk) * (1. - Yj(i));
        }
        LLj = log(arma::accu(pic % pis)) + LLj;
    }
    return LLj;
}

double rRUM_LLj(unsigned int N, unsigned int K, unsigned int J,
                const arma::vec &qj, const arma::vec &Yj, const arma::vec &Sj,
                const arma::vec &Gj, const arma::mat &alpha)
{
    double LLj = 0;
    for (unsigned int i = 0; i < N; i++) {
        double prod_ijk = 1.0;
        for (unsigned int k = 0; k < K; k++) {
            prod_ijk = (qj(k) * ((1. - Sj(k)) * alpha(i, k) +
                                 Gj(k) * (1. - alpha(i, k))) +
                        1. - qj(k)) *
                       prod_ijk;
        }
        LLj = log(prod_ijk * Yj(i) + (1. - prod_ijk) * (1. - Yj(i))) + LLj;
    }
    return LLj;
}

double nus_MH_lik(unsigned int N, unsigned int K, unsigned int nClass,
                  const arma::vec &Y_tm1, const arma::vec &Y_ast,
                  const arma::vec CLASS, const arma::vec &q_tm1,
                  const arma::vec &r_tm1, const arma::vec &r_ast, double pi_tm1,
                  double pi_ast)
{
    // create cubes of class by y by qj
    arma::cube lnYij_tm1(nClass, 2, 2); // cube for old row
    arma::cube lnYij_ast(nClass, 2, 2); // cube for candidate row
    for (unsigned int cc = 0; cc < nClass; cc++) {
        arma::vec ac = inv_bijectionvector(K, cc);
        double prod_ijk_tm1_ek = 1.0;
        double prod_ijk_tm1_1K = 1.0;
        double prod_ijk_ast_ek = 1.0;
        double prod_ijk_ast_1K = 1.0;
        // assuming qj_ast = 1_K to compute P(Y=1)
        for (unsigned int k = 0; k < K; k++) {
            prod_ijk_tm1_ek *= (q_tm1(k) * (1. - ac(k)) * r_tm1(k) + 1. -
                                q_tm1(k) * (1. - ac(k)));
            prod_ijk_tm1_1K *= ((1. - ac(k)) * r_tm1(k) + 1. - (1. - ac(k)));
            prod_ijk_ast_ek *= (q_tm1(k) * (1. - ac(k)) * r_ast(k) + 1. -
                                q_tm1(k) * (1. - ac(k)));
            prod_ijk_ast_1K *= ((1. - ac(k)) * r_ast(k) + 1. - (1. - ac(k)));
        }
        prod_ijk_tm1_ek = pi_tm1 * prod_ijk_tm1_ek;
        prod_ijk_tm1_1K = pi_tm1 * prod_ijk_tm1_1K;
        prod_ijk_ast_ek = pi_ast * prod_ijk_ast_ek;
        prod_ijk_ast_1K = pi_ast * prod_ijk_ast_1K;
        lnYij_tm1(cc, 0, 0) = log(1. - prod_ijk_tm1_ek);
        lnYij_tm1(cc, 1, 0) = log(prod_ijk_tm1_ek);
        lnYij_tm1(cc, 0, 1) = log(1. - prod_ijk_tm1_1K);
        lnYij_tm1(cc, 1, 1) = log(prod_ijk_tm1_1K);
        lnYij_ast(cc, 0, 0) = log(1. - prod_ijk_ast_ek);
        lnYij_ast(cc, 1, 0) = log(prod_ijk_ast_ek);
        lnYij_ast(cc, 0, 1) = log(1. - prod_ijk_ast_1K);
        lnYij_ast(cc, 1, 1) = log(prod_ijk_ast_1K);
    }
    double LL_old = 0.;
    double LL_new = 0.;
    for (unsigned int i = 0; i < N; i++) {
        double Y_tm1_i = Y_tm1(i);
        double Y_ast_i = Y_ast(i);
        double alpha_i = CLASS(i);
        LL_old +=
            lnYij_tm1(alpha_i, Y_tm1_i, 0) + lnYij_ast(alpha_i, Y_ast_i, 1);
        LL_new +=
            lnYij_tm1(alpha_i, Y_tm1_i, 1) + lnYij_ast(alpha_i, Y_ast_i, 0);
    }
    arma::vec ln_rMH(2);
    ln_rMH(0) = 0.;
    ln_rMH(1) = LL_new - LL_old;
    return min(ln_rMH);
}

double Beta(double a, double b)
{
    arma::vec avec = arma::linspace<arma::vec>(1, a - 1, a - 1);
    arma::vec bvec = arma::linspace<arma::vec>(b, a + b - 1, a);
    double lna = arma::accu(log(avec));
    double lnb = arma::accu(log(bvec));
    return exp(lna - lnb);
}

double ln_Beta(double a, double b)
{
    double flag = (1. * (a > 99)) * (1. * (b > 99));
    double ln_beta;
    if (flag == 0) {
        arma::vec avec = arma::linspace<arma::vec>(1, a - 1, a - 1);
        arma::vec bvec = arma::linspace<arma::vec>(b, a + b - 1, a);
        double lna = arma::accu(log(avec));
        double lnb = arma::accu(log(bvec));
        ln_beta = lna - lnb;
    } else {
        ln_beta = log(2. * arma::datum::pi) + (a - .5) * log(a) +
                  (b - .5) * log(b) - (a + b - .5) * log(a + b);
    }
    return ln_beta;
}

// spike/slab for sjk, gjk based on djk + nus for sampling 2 I_K rows + uses lik
// with Y

Rcpp::List parm_update11(unsigned int N, unsigned int J, unsigned int K,
                         unsigned int nClass, const arma::mat Y, arma::mat &Q,
                         arma::vec &CLASS, arma::mat &X, arma::mat &Smat,
                         arma::mat &Gmat, arma::vec &pis,
                         const arma::cube &match, const arma::mat &PmatYeq1,
                         const arma::vec &vv, arma::mat &nus, arma::mat &deltas,
                         double &omega)
{
    arma::cube Pcube(nClass, nClass, J);
    arma::cube abparms = arma::zeros<arma::cube>(J, K, 4);
    // Update X
    for (unsigned int j = 0; j < J; j++) {
        arma::rowvec qj = Q.row(j);
        double qj_biject = arma::conv_to<double>::from(qj * vv);
        arma::vec Pj_cl_Yeq1 = PmatYeq1.col(qj_biject);
        arma::rowvec gj = Gmat.row(j);
        arma::rowvec sj = Smat.row(j);
        arma::vec Yj = Y.col(j);
        arma::mat Pmatj = P_X_CL_Yeq0(K, nClass, qj.t(), gj.t(), sj.t());
        arma::rowvec denoms = arma::sum(Pmatj, 0);
        for (unsigned int i = 0; i < N; i++) {
            double Yij = Yj(i);
            double class_i = CLASS(i);
            if (Yij == 1) {
                X(i, j) = rgen::rmultinomial(Pj_cl_Yeq1);
            } else {
                arma::vec Pj_cl = Pmatj.col(class_i) / denoms(class_i);
                X(i, j) = rgen::rmultinomial(Pj_cl);
            }
        }
        Pcube.slice(j) = P_X_CL(K, nClass, qj.t(), gj.t(), sj.t()); // Pmatj;
    }
    // update CLASSES
    for (unsigned int i = 0; i < N; i++) {
        arma::vec pCL = arma::ones<arma::vec>(nClass);
        arma::rowvec Xi = X.row(i);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            for (unsigned int j = 0; j < J; j++) {
                double Xij = Xi(j);
                pCL(cc) *= Pcube(Xij, cc, j);
            }
        }
        arma::vec numerator = pCL % pis;
        arma::vec PS = numerator / arma::sum(numerator);
        double class_i = rgen::rmultinomial(PS);
        CLASS(i) = class_i;
        // update guess and slip full conditional beta parms
        for (unsigned int j = 0; j < J; j++) {
            arma::vec ab_ij = match.tube(Xi(j), class_i);
            for (unsigned int k = 0; k < K; k++) {
                abparms(j, k, ab_ij(k)) += 1.;
            }
        }
    }
    // update pis
    arma::uvec class_sum =
        arma::hist(CLASS, arma::linspace<arma::vec>(0, nClass - 1, nClass));
    arma::vec deltatilde = arma::conv_to<arma::vec>::from(class_sum) + 1.;
    pis = rgen::rdirichlet(deltatilde);
    // update Smat and Gmat, deltas
    double omega_new;
    arma::vec pistar = arma::zeros<arma::vec>(J);
    arma::mat rstar = arma::zeros<arma::mat>(J, K);
    double pg, ps, ug, us, gjk, sjk;
    // finding js to update via spike
    arma::vec no_nu_j = arma::ones<arma::vec>(J);
    for (unsigned int k = 0; k < K; k++) {
        for (unsigned int l = 0; l < 2; l++) {
            no_nu_j(nus(k, l)) = 0;
        }
    }
    for (unsigned int j = 0; j < J; j++) {
        double nu_flag = no_nu_j(j);
        double pistar_temp = 1.0;
        arma::rowvec qj = Q.row(j);
        for (unsigned int k = 0; k < K; k++) {
            double agk = abparms(j, k, 2);
            double bgk = abparms(j, k, 0);
            double bsk = abparms(j, k, 3);
            double ask = abparms(j, k, 1);
            // update deltajk
            double ln_pdjk;
            double qjk = qj(k);
            if (qjk == 1) {
                double m_sum = ask / (ask + bsk) + agk / (agk + bgk);
                double v_sum =
                    ask * bsk / (pow(ask + bsk, 2) * (ask + bsk + 1)) +
                    agk * bgk / (pow(agk + bgk, 2) * (agk + bgk + 1));
                double ln_pnorm = R::pnorm(0.0, m_sum - 1., sqrt(v_sum), 1,
                                           1); // last '1' denotes log p
                double ln_pX_d_eq_1 = ln_pnorm + ln_Beta(ask + 1., bsk + 1.) +
                                      ln_Beta(agk + 1., bgk + 1.);
                double ln_pX_d_eq_0 = ln_Beta(bsk + agk + 1, ask + bgk + 1);
                ln_pdjk =
                    ln_pX_d_eq_0 + log(1. - omega) - log(omega) - ln_pX_d_eq_1;
            } else {
                ln_pdjk = log(1. - omega) - log(omega);
            }
            double ud = R::runif(0.0, 1.0);
            double djk = (log(1. - ud) - log(ud) > ln_pdjk) * nu_flag +
                         (1. - nu_flag) * qjk;
            deltas(j, k) = djk;
            // deltajk=1 => g<1-s
            ug = R::runif(0.0, 1.0);
            if (djk == 1) {
                us = R::runif(0.0, 1.0);
                if (qjk == 0) {
                    gjk = 1. - sqrt(1. - ug);
                    sjk = us * (1. - gjk);
                } else {
                    // draw g conditoned upon s_t-1
                    pg = R::pbeta(1.0 - Smat(j, k), agk + 1.0, bgk + 1.0, 1, 0);
                    gjk = R::qbeta(ug * pg, agk + 1.0, bgk + 1.0, 1, 0);
                    // draw s conditoned upon g
                    ps = R::pbeta(1.0 - gjk, ask + 1.0, bsk + 1.0, 1, 0);
                    sjk = R::qbeta(us * ps, ask + 1.0, bsk + 1.0, 1, 0);
                }
            }
            // deltajk=0 => g=1-s
            else {
                if (qjk == 0) {
                    gjk = ug * nu_flag + (1. - nu_flag) * .5;
                    sjk = 1. - gjk;
                } else {
                    gjk = R::qbeta(ug, bsk + agk + 1.0, ask + bgk + 1.0, 1, 0);
                    sjk = 1. - gjk;
                }
            }
            Gmat(j, k) = gjk;
            Smat(j, k) = sjk;
            rstar(j, k) = gjk / (1.0 - sjk); // compute rstarjk
            if (qjk == 1.) {
                pistar_temp = (1.0 - sjk) * pistar_temp; // compute pistarj
            }
            // update omega
            double sum_deltas =
                arma::accu(deltas); //-2.*K;//subtracting those from nu rows
            double uw = R::runif(0.0, 1.0);
            omega_new =
                R::qbeta(uw, sum_deltas + 1.0, J * K - sum_deltas + 1.0, 1, 0);
        }
        pistar(j) = pistar_temp;
    }
    // Update Q
    arma::vec one_K = arma::ones<arma::vec>(K);
    arma::vec one_to_K = arma::linspace(0, K - 1, K);
    for (unsigned int k = 0; k < K; k++) {
        arma::uvec no_k = find(one_to_K != k);
        arma::vec ek = inv_bijectionvector(K, vv(k));
        // looping over l for nuks
        for (unsigned int l = 0; l < 2; l++) {
            arma::vec valid_k = arma::ones<arma::vec>(J);
            valid_k(nus(k, 0)) = 0;
            valid_k(nus(k, 1)) = 0;
            for (unsigned int kk = 0; kk < no_k.n_elem; kk++) {
                for (unsigned int s = 0; s < 2; s++) {
                    double j = nus(no_k(kk), s);
                    valid_k(j) = 0;
                }
            }
            arma::uvec temp_ones = find(valid_k == 1);
            unsigned int J_temp = temp_ones.n_elem;
            unsigned int j = arma::conv_to<unsigned int>::from(
                arma::randi(1, arma::distr_param(0, J_temp - 1)));
            double j_ast = temp_ones(j);
            double j_nukl = nus(k, l);
            double ln_u = log(R::runif(0.0, 1.0));
            double min_r =
                nus_MH_lik(N, K, nClass, Y.col(j_nukl), Y.col(j_ast), CLASS, ek,
                           (rstar.row(j_nukl)).t(), (rstar.row(j_ast)).t(),
                           pistar(j_nukl), pistar(j_ast));
            double select = (min_r > ln_u);
            if (select == 1.) {
                double nu_new = j_ast * select + (1. - select) * j_nukl;
                nus(k, l) = nu_new;
                Q.row(j_ast) = ek.t();
                Q.row(j_nukl) = one_K.t();
            }
        }
    }

    return Rcpp::List::create(Rcpp::Named("pistar", pistar),
                              Rcpp::Named("rstar", rstar),
                              Rcpp::Named("omega_new", omega_new));
}

// spike/slab for sjk, gjk based on djk + nus for sampling 2 I_K rows + uses lik
// with Y

Rcpp::List parm_update12(unsigned int N, unsigned int J, unsigned int K,
                         unsigned int nClass, const arma::mat &Y, arma::mat &Q,
                         arma::vec &CLASS, arma::mat &X, arma::mat &Smat,
                         arma::mat &Gmat, arma::vec &pis,
                         const arma::cube &match, const arma::mat &PmatYeq1,
                         const arma::vec &vv, arma::mat &nus, arma::mat &deltas,
                         double &omega)
{
    arma::cube Pcube(nClass, nClass, J);
    arma::cube abparms = arma::zeros<arma::cube>(J, K, 4);
    // Update X
    for (unsigned int j = 0; j < J; j++) {
        arma::rowvec qj = Q.row(j);
        double qj_biject = arma::conv_to<double>::from(qj * vv);
        arma::vec Pj_cl_Yeq1 = PmatYeq1.col(qj_biject);
        arma::rowvec gj = Gmat.row(j);
        arma::rowvec sj = Smat.row(j);
        arma::vec Yj = Y.col(j);
        arma::mat Pmatj = P_X_CL_Yeq0(K, nClass, qj.t(), gj.t(), sj.t());
        arma::rowvec denoms = arma::sum(Pmatj, 0);
        for (unsigned int i = 0; i < N; i++) {
            double Yij = Yj(i);
            double class_i = CLASS(i);
            if (Yij == 1) {
                X(i, j) = rgen::rmultinomial(Pj_cl_Yeq1);
            } else {
                arma::vec Pj_cl = Pmatj.col(class_i) / denoms(class_i);
                X(i, j) = rgen::rmultinomial(Pj_cl);
            }
        }
        Pcube.slice(j) = P_X_CL(K, nClass, qj.t(), gj.t(), sj.t());
    }
    // update CLASSES
    for (unsigned int i = 0; i < N; i++) {
        arma::vec pCL = arma::ones<arma::vec>(nClass);
        arma::rowvec Xi = X.row(i);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            for (unsigned int j = 0; j < J; j++) {
                double Xij = Xi(j);
                pCL(cc) *= Pcube(Xij, cc, j);
            }
        }
        arma::vec numerator = pCL % pis;
        arma::vec PS = numerator / arma::sum(numerator);
        double class_i = rgen::rmultinomial(PS);
        CLASS(i) = class_i;
        // update guess and slip full conditional beta parms
        for (unsigned int j = 0; j < J; j++) {
            arma::vec ab_ij = match.tube(Xi(j), class_i);
            for (unsigned int k = 0; k < K; k++) {
                abparms(j, k, ab_ij(k)) += 1.;
            }
        }
    }
    // update pis
    arma::uvec class_sum =
        arma::hist(CLASS, arma::linspace<arma::vec>(0, nClass - 1, nClass));
    arma::vec deltatilde = arma::conv_to<arma::vec>::from(class_sum) + 1.;
    pis = rgen::rdirichlet(deltatilde);
    // update Smat and Gmat, deltas
    double omega_new;
    arma::vec pistar = arma::zeros<arma::vec>(J);
    arma::mat rstar = arma::zeros<arma::mat>(J, K);
    double pg, ps, ug, us, gjk, sjk;
    // finding js to update via spike
    arma::vec no_nu_j = arma::ones<arma::vec>(J);
    for (unsigned int k = 0; k < K; k++) {
        for (unsigned int l = 0; l < 2; l++) {
            no_nu_j(nus(k, l)) = 0;
        }
    }
    for (unsigned int j = 0; j < J; j++) {
        double nu_flag = no_nu_j(j);
        double pistar_temp = 1.0;
        arma::rowvec qj = Q.row(j);
        for (unsigned int k = 0; k < K; k++) {
            double agk = abparms(j, k, 2);
            double bgk = abparms(j, k, 0);
            double bsk = abparms(j, k, 3);
            double ask = abparms(j, k, 1);
            // update deltajk
            double ln_pdjk, djk;
            double qjk = qj(k);
            if (nu_flag == 1) {
                double m_sum = ask / (ask + bsk) + agk / (agk + bgk);
                double v_sum =
                    ask * bsk / (pow(ask + bsk, 2) * (ask + bsk + 1)) +
                    agk * bgk / (pow(agk + bgk, 2) * (agk + bgk + 1));
                double ln_pnorm = R::pnorm(0.0, m_sum - 1., sqrt(v_sum), 1,
                                           1); // last '1' denotes log p
                double ln_pX_d_eq_1 = ln_pnorm + ln_Beta(ask + 1., bsk + 1.) +
                                      ln_Beta(agk + 1., bgk + 1.);
                double ln_pX_d_eq_0 = ln_Beta(bsk + agk + 1, ask + bgk + 1);
                ln_pdjk =
                    ln_pX_d_eq_0 + log(1. - omega) - log(omega) - ln_pX_d_eq_1;
                double ud = R::runif(0.0, 1.0);
                djk = (log(1. - ud) - log(ud) > ln_pdjk);
                deltas(j, k) = djk;
                ug = R::runif(0.0, 1.0);
                if (djk == 1) {
                    us = R::runif(0.0, 1.0);
                    // draw g conditoned upon s_t-1
                    pg = R::pbeta(1.0 - Smat(j, k), agk + 1.0, bgk + 1.0, 1, 0);
                    gjk = R::qbeta(ug * pg, agk + 1.0, bgk + 1.0, 1, 0);
                    // draw s conditoned upon g
                    ps = R::pbeta(1.0 - gjk, ask + 1.0, bsk + 1.0, 1, 0);
                    sjk = R::qbeta(us * ps, ask + 1.0, bsk + 1.0, 1, 0);
                }
                // deltajk=0 => g=1-s
                else {
                    gjk = R::qbeta(ug, bsk + agk + 1.0, ask + bgk + 1.0, 1, 0);
                    sjk = 1. - gjk;
                }
            } else {
                djk = qjk;
                deltas(j, k) = djk;
                if (qjk == 0) {
                    gjk = .5;
                    sjk = 1. - gjk;
                } else {
                    ug = R::runif(0.0, 1.0);
                    us = R::runif(0.0, 1.0);
                    pg = R::pbeta(1.0 - Smat(j, k), agk + 1.0, bgk + 1.0, 1,
                                  0); // draw g conditoned upon s_t-1
                    gjk = R::qbeta(ug * pg, agk + 1.0, bgk + 1.0, 1, 0);
                    ps = R::pbeta(1.0 - gjk, ask + 1.0, bsk + 1.0, 1,
                                  0); // draw s conditoned upon g
                    sjk = R::qbeta(us * ps, ask + 1.0, bsk + 1.0, 1, 0);
                }
            }
            Gmat(j, k) = gjk;
            Smat(j, k) = sjk;
            rstar(j, k) = gjk / (1.0 - sjk); // compute rstarjk
            if (qjk == 1.) {
                pistar_temp = (1.0 - sjk) * pistar_temp; // compute pistarj
            }
        }
        pistar(j) = pistar_temp;
    }
    // update omega
    double sum_deltas =
        arma::accu(deltas) - 2. * K; // subtracting ones from nu rows
    double uw = R::runif(0.0, 1.0);
    omega_new = R::qbeta(uw, sum_deltas + 1.0,
                         (J - 2 * K) * K - sum_deltas + 1.0, 1, 0);
    // Update Q
    arma::vec one_K = arma::ones<arma::vec>(K);
    arma::vec one_to_K = arma::linspace(0, K - 1, K);
    for (unsigned int k = 0; k < K; k++) {
        arma::uvec no_k = find(one_to_K != k);
        arma::vec ek = inv_bijectionvector(K, vv(k));
        // looping over l for nuks
        for (unsigned int l = 0; l < 2; l++) {
            arma::vec valid_k = arma::ones<arma::vec>(J);
            valid_k(nus(k, 0)) = 0;
            valid_k(nus(k, 1)) = 0;
            for (unsigned int kk = 0; kk < no_k.n_elem; kk++) {
                for (unsigned int s = 0; s < 2; s++) {
                    double j = nus(no_k(kk), s);
                    valid_k(j) = 0;
                }
            }
            arma::uvec temp_ones = find(valid_k == 1);
            unsigned int J_temp = temp_ones.n_elem;
            unsigned int j = arma::conv_to<unsigned int>::from(
                arma::randi(1, arma::distr_param(0, J_temp - 1)));
            double j_ast = temp_ones(j);
            double j_nukl = nus(k, l);
            double ln_u = log(R::runif(0.0, 1.0));
            double min_r =
                nus_MH_lik(N, K, nClass, Y.col(j_nukl), Y.col(j_ast), CLASS, ek,
                           (rstar.row(j_nukl)).t(), (rstar.row(j_ast)).t(),
                           pistar(j_nukl), pistar(j_ast));
            double select = (min_r > ln_u);
            if (select == 1.) {
                double nu_new = j_ast * select + (1. - select) * j_nukl;
                nus(k, l) = nu_new;
                Q.row(j_ast) = ek.t();
                Q.row(j_nukl) = one_K.t();
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("pistar", pistar),
                              Rcpp::Named("rstar", rstar),
                              Rcpp::Named("omega_new", omega_new));
}

// version that sets Q = 1, spike slab prior, swaps 2I_K rows with MH steps

// [[Rcpp::export]]
Rcpp::List errum_Gibbs_Q(const arma::mat &Y, unsigned int K, unsigned int burnin,
                         unsigned int chain_length = 10000)
{
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int nClass = pow(2, K);

    // Total Chain length with burn
    // Prevents overflow
    unsigned int iter_total = chain_length + burnin;

    unsigned int tmburn;
    arma::vec vv = bijectionvector(K);

    // Saving output
    arma::mat PISTAR(J, chain_length);
    arma::cube RSTAR(J, K, chain_length);
    arma::mat PIs(nClass, chain_length);
    arma::cube QS(J, K, chain_length);
    arma::mat m_Delta = arma::zeros<arma::mat>(J, K);
    arma::mat Delta_biject(J, chain_length);
    arma::cube NUS(K, 2, chain_length);

    // need to initialize, alphas, X,ss, gs,pis
    arma::vec CLASS =
        arma::randi<arma::vec>(N, arma::distr_param(0, nClass - 1));
    arma::mat X = arma::zeros<arma::mat>(N, J);
    arma::mat ss = arma::randu<arma::mat>(J, K);
    arma::mat gs =
        (arma::ones<arma::mat>(J, K) - ss) % arma::randu<arma::mat>(J, K);
    arma::vec delta0 = arma::ones<arma::vec>(nClass);
    arma::vec pis = rgen::rdirichlet(delta0);
    arma::mat Q = arma::ones<arma::mat>(J, K);
    arma::vec one_to_J = arma::linspace(0, J - 1, J);
    arma::vec one_to_J_shuffled = arma::shuffle(one_to_J);
    arma::mat nus(K, 2);
    double ind = 0;
    for (unsigned int m1 = 0; m1 < 2; m1++) {
        for (unsigned int m2 = 0; m2 < K; m2++) {
            nus(m2, m1) = one_to_J_shuffled(ind);
            Q.row(one_to_J_shuffled(ind)) =
                (inv_bijectionvector(K, vv(m2))).t();
            ind += 1;
        }
    }
    arma::cube match = X_CL_match(K, nClass);
    arma::mat PmatYeq1 = P_X_CL_Yeq1(K, nClass);
    arma::mat deltas = arma::randu<arma::mat>(J, K); // K>1 is assumed
    deltas.elem(find(deltas > 0.5)).ones();
    deltas.elem(find(deltas <= 0.5)).zeros();
    double omega = R::runif(0, 1);

    arma::mat M1(chain_length, J);
    arma::cube M2(J, J, chain_length);

    // Start Markov chain
    for (unsigned int t = 0; t < iter_total; t++) {
        Rcpp::List output =
            parm_update12(N, J, K, nClass, Y, Q, CLASS, X, ss, gs, pis, match,
                          PmatYeq1, vv, nus, deltas, omega);
        omega = Rcpp::as<double>(output[2]);
        if (t > burnin - 1) {
            tmburn = t - burnin;
            arma::vec pistar = Rcpp::as<arma::vec>(output[0]);
            PISTAR.col(tmburn) = pistar;
            arma::mat rstar = Rcpp::as<arma::mat>(output[1]);
            RSTAR.slice(tmburn) = rstar;
            PIs.col(tmburn) = pis;
            QS.slice(tmburn) = Q;
            m_Delta = (deltas + m_Delta * tmburn) / (tmburn + 1.);
            Delta_biject.col(tmburn) = deltas * vv;
            // refactor
            arma::mat Yt = simrRUM_errum(N, J, K, Q, rstar, pistar, CLASS);
            M1.row(tmburn) = arma::mean(Yt);
            M2.slice(tmburn) = OddsRatio(N, J, Yt);
            NUS.slice(tmburn) = nus;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("PISTAR", PISTAR), Rcpp::Named("RSTAR", RSTAR),
        Rcpp::Named("PIs", PIs), Rcpp::Named("QS", QS),
        Rcpp::Named("m_Delta", m_Delta),
        Rcpp::Named("Delta_biject", Delta_biject), Rcpp::Named("M2", M2),
        Rcpp::Named("M1", M1), Rcpp::Named("NUS", NUS));
}
