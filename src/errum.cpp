#include "ecdm.h"
#include <rgen.h>

//' Simulate rRUM Data
//'
//' Generate data that falls underneath a rRUM modeling paradigm.
//'
//' @param N        Number of Subjects
//' @param J        Number of Assessment Items
//' @param K        Number of Attributes
//' @param rstar    Matrix of penalities for missing an attribute on an item
//'                 with dimensions \eqn{J \times K}{J x K}.
//' @param pistar   Vector of probabilities associated with whether
//'                 \eqn{Y_{ij} = 1} if subject \eqn{i} has all the required
//'                 attributes having length \eqn{J}.
//' @param Q        Matrix of attributes required for each assessment
//'                 item with dimensions \eqn{J \times K}{J x K}.
//' @param alpha    Matrix with dimensions \eqn{N \times K}{N x K} of 1s and 0s.
//' @export
// [[Rcpp::export]]
arma::mat simrRUM(unsigned int N, unsigned int J, unsigned int K,
                  const arma::mat& Q, const arma::mat& rstar,
                  const arma::vec& pistar, const arma::mat& alpha){

    arma::vec k_index = arma::linspace(0,K-1,K);
    double kj;
    double aik;
    arma::vec pmu(J);
    arma::mat Y=arma::zeros<arma::mat>(N,J);
    for(unsigned int i=0; i<N;i++){
        arma::vec Yi=arma::zeros<arma::vec>(J);
        arma::vec pi=arma::ones<arma::vec>(J);
        arma::vec ui = arma::randu<arma::vec>(J);
        for(unsigned int j=0;j<J;j++){
            arma::uvec task_ij = find(Q.row(j) == 1);

            for(unsigned int k = 0;k<task_ij.n_elem ;k++){
                kj = task_ij(k);
                aik = alpha(i,kj);
                pi(j) = ((rstar(j,kj)*(1.0-aik)+1.0*aik)*Q(j,kj)+1.0*(1.0-Q(j,kj)))*pi(j);
            }
            pi(j) = pistar(j)*pi(j);
        }
        pmu = pi - ui;
        Yi(arma::find(pmu > 0)).fill(1);
        Y.row(i) = Yi.t();
    }
    return Y;
}


//' Parameter Update Routine for rRUM Gibbs Sampler
//'
//' @param N        Number of Subjects
//' @param J        Number of Assessment Items
//' @param K        Number of Attributes
//' @param C        Number of Classes (\eqn{2^K})
//' @param Y        Matrix of binary responses to assessment items with
//'                 dimensions \eqn{N \times J}{N x J}.
//' @param Q        Matrix of attributes required for each assessment
//'                 item with dimensions \eqn{J \times K}{J x K}.
//' @param alpha    Matrix with dimensions \eqn{N \times K}{N x K} of 1s and 0s.
//' @param X        Cube with dimensions \eqn{N \times J \times K}{N x J x K}.
//' @param Smat     Slipping parameters
//' @param Gmat     Guessing Parameters
//' @param pi       Vector of PIs
//' @param vv       Bijection vector of attributes
//' @param delta0   Vector of length \eqn{2^K} that corresponds to the
//'                 prior values for betas and dirichlet distributions of the
//'                 classes.
//' @export
// [[Rcpp::export]]
Rcpp::List parm_update(unsigned int N, unsigned int J, unsigned int K,
                       unsigned int C,
                       const arma::mat Y, const arma::mat& Q,
                       arma::mat& alpha, arma::cube& X, arma::mat& Smat,
                       arma::mat& Gmat, arma::vec& pi,
                       const arma::vec vv, const arma::vec& delta0){

    arma::vec k_index = arma::linspace(0,K-1,K);
    double kj,prodXijk,pi_ijk,aik,u,compare;
    double pi_ik,aik_nmrtr_k,aik_dnmntr_k,c_aik_1,c_aik_0;
    arma::vec aik_nmrtr(K);
    arma::vec aik_dnmntr(K);

    //Update X cube
    for(unsigned int i=0;i<N;i++){

        arma::vec ai = (alpha.row(i)).t();
        arma::vec Yi = (Y.row(i)).t();
        arma::vec ui = arma::randu<arma::vec>(K);
        aik_nmrtr    = arma::ones<arma::vec>(K);
        aik_dnmntr   = arma::ones<arma::vec>(K);

        for(unsigned int j=0;j<J;j++){

            double Yij = Yi(j);
            arma::vec Xij = X.tube(i,j);
            arma::uvec task_ij = find(Q.row(j) == 1);

            for(unsigned int k = 0;k<task_ij.n_elem ;k++){
                kj = task_ij(k);
                aik = alpha(i,kj);
                Xij(kj) = 1;
                prodXijk = prod(Xij(task_ij));
                u = R::runif(0.0,1.0);
                pi_ijk = (1.0-prodXijk)*(aik*(1.0-Smat(j,kj)) + (1.0-aik)*Gmat(j,kj) );
                compare=(pi_ijk>u);
                Xij(kj)=(1.0-Yij)*compare + Yij;

                aik_nmrtr(kj) = ( Xij(kj)*(1.0-Smat(j,kj)) + (1.0-Xij(kj))*Smat(j,kj) )*aik_nmrtr(kj);
                aik_dnmntr(kj) = ( Xij(kj)*Gmat(j,kj) + (1.0-Xij(kj))*(1.0-Gmat(j,kj)) )*aik_dnmntr(kj);
            }
            X.tube(i,j) = Xij;
        }

        //Update alpha_ik
        for(unsigned int k=0;k<K;k++){
            ai(k) = 1.0;
            c_aik_1 = (arma::as_scalar( ai.t()*vv ));
            ai(k) = 0.0;
            c_aik_0 = (arma::as_scalar( ai.t()*vv ));

            aik_nmrtr_k = aik_nmrtr(k)*pi(c_aik_1);
            aik_dnmntr_k = aik_dnmntr(k)*pi(c_aik_0);
            pi_ik = aik_nmrtr_k/(aik_nmrtr_k + aik_dnmntr_k);
            ai(k) = 1.0*(pi_ik > ui(k));
        }
        alpha.row(i) = ai.t();
    }

    //update pi
    arma::vec a_bijection = alpha * vv;
    arma::uvec deltatilde = arma::hist( a_bijection,arma::linspace<arma::vec>(0,C-1,C) );
    pi = rgen::rdirichlet(deltatilde+delta0);

    //update Smat and Gmat
    arma::vec pistar = arma::zeros<arma::vec>(J);
    arma::mat rstar = arma::zeros<arma::mat>(J,K);
    double pg,ps,ug,us,gjk,sjk;

    for(unsigned int j=0;j<J;j++){
        arma::uvec task_ij = find(Q.row(j) == 1);
        arma::mat Xj = X.tube(0,j,N-1,j);
        double pistar_temp =1.0;

        for(unsigned int k = 0;k<task_ij.n_elem ;k++){
            kj = task_ij(k);
            arma::vec Xjk = Xj.col(kj);
            arma::vec ak = alpha.col(kj);

            double Sumalphak =  (arma::as_scalar(ak.t() * ak));
            double SumXjk = (arma::as_scalar(Xjk.t() * Xjk));
            double SumXjkalphak = (arma::as_scalar(Xjk.t() * ak));
            double bsk = SumXjkalphak ;
            double ask = Sumalphak - SumXjkalphak ;
            double agk = SumXjk - SumXjkalphak ;
            double bgk = N - SumXjk - Sumalphak + SumXjkalphak ;
            ug = R::runif(0.0,1.0);
            us = R::runif(0.0,1.0);

            //draw g conditoned upon s_t-1
            pg = R::pbeta(1.0-Smat(j,kj),agk+1.0,bgk+1.0,1,0);
            gjk = R::qbeta(ug*pg,agk+1.0,bgk+1.0,1,0);
            //draw s conditoned upon g
            ps = R::pbeta(1.0-gjk,ask+1.0,bsk+1.0,1,0);
            sjk = R::qbeta(us*ps,ask+1.0,bsk+1.0,1,0);

            Gmat(j,kj) = gjk;
            Smat(j,kj) = sjk;

            rstar(j,kj) = gjk/(1.0 - sjk);//compute rstarjk
            pistar_temp = (1.0-sjk)*pistar_temp;//compute pistarj
        }
        pistar(j) = pistar_temp;
    }

    return Rcpp::List::create(Rcpp::Named("pistar",pistar),
                              Rcpp::Named("rstar",rstar)
    );
}


//' rRUM Gibbs-based Estimation
//'
//' Performs an estimation of the rRUM model using the Metropolis-Hastings
//' derivation.
//'
//' @param Y            Matrix of binary responses to assessment items with
//'                     dimensions \eqn{N \times J}{N x J}.
//' @param Q            Matrix of attributes required for each assessment
//'                     item with dimensions \eqn{J \times K}{J x K}.
//' @param chain_length Number of iterations
//' @export
// [[Rcpp::export]]
Rcpp::List rRUM_Gibbs(const arma::mat& Y, const arma::mat& Q,
                      unsigned int chain_length = 10000){
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int K = Q.n_cols;
    unsigned int C = pow(2,K);

    arma::vec vv = bijectionvector(K);

    //Prior values for betas and Dirichlet distribution
    arma::vec delta0 = arma::ones<arma::vec>(C);

    //Savinging output
    arma::mat PISTAR(J,chain_length);
    arma::cube RSTAR(J,K,chain_length);
    arma::mat PIs(C,chain_length);
    //arma::mat ALPHAS(N,chain_length);

    //need to initialize, alphas, X,ss, gs,pis
    arma::mat alpha = arma::randu<arma::mat>(N,K); //K>1 is assumed
    alpha.elem( find(alpha > 0.5) ).ones();
    alpha.elem( find(alpha <= 0.5) ).zeros();
    arma::cube X = arma::zeros<arma::cube>(N,J,K);
    arma::mat ss = arma::randu<arma::mat>(J,K);
    arma::mat gs = (arma::ones<arma::mat>(J,K) - ss)%arma::randu<arma::mat>(J,K);
    arma::vec pis = rgen::rdirichlet(delta0);

    //Start Markov chain
    for(unsigned int t = 0; t < chain_length; t++){
        //updata X,alpha,pi,s,g,pistar,rstar
        Rcpp::List output = parm_update(N,J,K,C,Y,Q,alpha,X,ss,gs,pis,vv,delta0);

        //update value for pis. alphas are updated via pointer. save classes and PIs
        PISTAR.col(t)  = Rcpp::as<arma::vec>(output[0]);
        RSTAR.slice(t) = Rcpp::as<arma::mat>(output[1]);
        PIs.col(t)     = pis;
        //ALPHAS.col(t)  = alpha*vv;
    }
    return Rcpp::List::create(Rcpp::Named("PISTAR",PISTAR),
                              Rcpp::Named("RSTAR",RSTAR),
                              Rcpp::Named("PIs",PIs)//,
                              //Rcpp::Named("ALPHAS",ALPHAS)
    );
}

//' Update Parameters in Metropolis-Hastings Approach
//'
//' Helper functions that handles the process of updating the parameters
//' under the derving MH approach.
//'
//' @param N        Number of Subjects
//' @param J        Number of Assessment Items
//' @param K        Number of Attributes
//' @param C        Number of Classes (\eqn{2^K})
//' @param Y        Matrix of binary responses to assessment items with
//'                 dimensions \eqn{N \times J}{N x J}.
//' @param Q        Matrix of attributes required for each assessment
//'                 item with dimensions \eqn{J \times K}{J x K}.
//' @param alpha    Matrix with dimensions \eqn{N \times K}{N x K} of 1s and 0s.
//' @param pistar   Vector of probabilities associated with whether
//'                 \eqn{Y_{ij} = 1} if subject \eqn{i} has all the required
//'                 attributes having length \eqn{J}.
//' @param rstar    Matrix with dimensions \eqn{J \times K}{J x K}.
//' @param pi       Vector of length \eqn{2^K} containing the
//'                 latent class probabilities that \eqn{Y_{ij} = 1}.
//' @param vv       Bijection vector of length \eqn{k} for attributes.
//' @param delta0   Vector of length \eqn{2^K} that corresponds to the
//'                 prior values for betas and dirichlet distributions of the
//'                 classes.
//' @param delta    TBA?
//' @param Amat     TBA?
//' @export
// [[Rcpp::export]]
Rcpp::List parm_update_MH(unsigned int N, unsigned int J,
                          unsigned int K, unsigned int C,
                          const arma::mat Y, const arma::mat& Q,
                          arma::mat& alpha, arma::vec& pistar, arma::mat& rstar,
                          arma::vec& pi, const arma::vec vv, const arma::vec& delta0,
                          double delta, const arma::mat& Amat){

    arma::vec k_index = arma::linspace(0,K-1,K);
    double kj;
    // double pi_ik,aik_nmrtr_k,aik_dnmntr_k,c_aik_1,c_aik_0;
    arma::vec aik_nmrtr(K);
    arma::vec aik_dnmntr(K);
    arma::mat Alpha = arma::zeros<arma::mat>(J,K);
    arma::mat Ac = arma::zeros<arma::mat>(J,K);
    arma::mat pistar_num = arma::ones<arma::mat>(N,J);
    arma::mat pistar_denom = arma::ones<arma::mat>(N,J);

    arma::mat PENALTY = arma::zeros<arma::mat>(N,J);
    arma::mat PISTARrep = arma::zeros<arma::mat>(N,J);
    arma::mat like = arma::zeros<arma::mat>(N,C);
    arma::vec Yi = arma::zeros<arma::vec>(J);

    arma::vec pistar_new = arma::zeros<arma::vec>(J);
    arma::mat rstar_new = arma::zeros<arma::mat>(J,K);
    double a,b,c,d;

    for(unsigned int j=0;j<J;j++){
        arma::uvec task_ij = find(Q.row(j) == 1);
        a = pistar(j) - delta;
        b = pistar(j) + delta;
        pistar_new(j) = R::runif(a, b);
        while((pistar_new(j) <= 0) || (pistar_new(j) >= 1)) pistar_new(j) = R::runif(a, b);
        for(unsigned int k = 0;k<task_ij.n_elem ;k++){
            kj = task_ij(k);
            c = rstar(j,kj) - delta;
            d = rstar(j,kj) + delta;
            rstar_new(j,kj) = R::runif(c,d);
            while((rstar_new(j,kj) <= 0) || (rstar_new(j,kj) >= 1)) rstar_new(j,kj) = R::runif(c,d);
        }
    }

    for(unsigned int c=0; c<C;c++){
        Ac= repmat(Amat.row(c), J, 1);
        PISTARrep = repmat(pistar.t(), N, 1);
        PENALTY = repmat(prod((rstar % (1-Ac) % Q) + (Ac % Q) + (1 - Q),1).t(), N, 1);
        like.col(c) = prod(PISTARrep % PENALTY % Y + (1 - PISTARrep % PENALTY) % (1 - Y), 1) * pi(c);
    }

    for(unsigned int i=0;i<N;i++){

        //Update alpha_ik
        alpha.row(i) = Amat.row(rgen::rmultinomial(like.row(i)/sum(like.row(i))));

        //calculate likelihoods for sampling pistar & rstar
        Yi = Y.row(i).t();
        Alpha = repmat(alpha.row(i), J, 1);
        pistar_denom.row(i) = ((pistar % prod((rstar % (1-Alpha) % Q) + (Q % Alpha) + (1 - Q),1) % Yi) +
            ((1-(pistar % prod((rstar % (1-Alpha) % Q) + (Q % Alpha) + (1 - Q),1))) % (1-Yi))).t();
pistar_num.row(i) = ((pistar_new % prod((rstar_new % (1-Alpha) % Q) + (Q % Alpha) + (1 - Q),1) % Yi) +
    ((1-(pistar_new % prod((rstar_new % (1-Alpha) % Q) + (Q % Alpha) + (1 - Q),1))) % (1-Yi))).t();
// double asdf = prod(pistar_num.row(i));

    }

    //update pi
    arma::vec a_bijection = alpha * vv;
    arma::uvec deltatilde = arma::hist( a_bijection,arma::linspace<arma::vec>(0,C-1,C) );
    pi = rgen::rdirichlet(deltatilde+delta0);

    //update pistar and rstar
    double phi = prod(vectorise(pistar_num/pistar_denom));
    double ar_pistar = 0;
    double ar_rstar = 0;
    double l = R::runif(0,1);
    if(l < phi){
        pistar = pistar_new;
        ar_pistar = 1;
    }
    l = R::runif(0,1);
    if(l < phi){
        rstar = rstar_new;
        ar_rstar = 1;
    }



    return Rcpp::List::create(Rcpp::Named("pistar",pistar),
                              Rcpp::Named("rstar",rstar),
                              Rcpp::Named("ar_rstar",ar_rstar),
                              Rcpp::Named("ar_pistar",ar_pistar)
    );
}


//' rRUM Metropolis-Hastings-based Estimation
//'
//' Performs an estimation of the rRUM model using the Metropolis-Hastings
//' derivation.
//'
//' @param Y            Matrix of binary responses to assessment items with
//'                     dimensions \eqn{N \times J}{N x J}.
//' @param Q            Matrix of attributes required for each assessment
//'                     item with dimensions \eqn{J \times K}{J x K}.
//' @param Amat         TBA?
//' @param delta        TBA?
//' @param chain_length Number of iterations
//' @return
//' A `list` containing:
//' - **PISTAR**: Estimates of the latent probability class.
//' - **RSTAR**: Estimates of the penalty parameter for missing a required \eqn{k}.
//' - **PIs**: Latent Class Probabilities with length \eqn{K}
//' - **ACCEPTREJECT_PISTAR**: Vector of 1's or 0's indicating whether the
//'                            values for pistar was accepted or rejected.
//' - **ACCEPTREJECT_RSTAR**: Vector of 1's or 0's indicating whether the
//'                            values for rstar was accepted or rejected.
//' @export
// [[Rcpp::export]]
Rcpp::List rRUM_MH(const arma::mat& Y, const arma::mat& Q,
                   arma::mat& Amat,
                   double delta, unsigned int chain_length = 10000){
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int K = Q.n_cols;
    unsigned int C = pow(2,K);

    arma::vec vv = bijectionvector(K);

    // Prior values for betas and Dirichlet distribution
    arma::vec delta0 = arma::ones<arma::vec>(C);

    //Savinging output
    arma::mat PISTAR(J,chain_length);
    arma::cube RSTAR(J,K,chain_length);
    arma::mat PIs(C,chain_length);
    arma::vec ACCEPTREJECT_PISTAR(chain_length);
    arma::vec ACCEPTREJECT_RSTAR(chain_length);
    //arma::mat ALPHAS(N,chain_length);

    //need to initialize, alphas, X,pistar, rstar,pis
    arma::mat alpha = arma::randu<arma::mat>(N,K); //K>1 is assumed
    alpha.elem( find(alpha > 0.5) ).ones();
    alpha.elem( find(alpha <= 0.5) ).zeros();
    arma::vec pistars = arma::randu<arma::vec>(J);
    arma::mat rstars = arma::randu<arma::mat>(J,K);
    arma::vec pis = rgen::rdirichlet(delta0);

    //Start Markov chain
    for(unsigned int t = 0; t < chain_length; t++){
        //updata X,alpha,pi,s,g,pistar,rstar
        Rcpp::List output = parm_update_MH(N, J, K, C, Y,
                                           Q, alpha, pistars, rstars,
                                           pis, vv, delta0, delta, Amat);

        //update value for pis. alphas are updated via pointer. save classes and PIs
        PISTAR.col(t)  = Rcpp::as<arma::vec>(output[0]);
        RSTAR.slice(t) = Rcpp::as<arma::mat>(output[1]);
        ACCEPTREJECT_RSTAR(t) = Rcpp::as<double>(output[2]);
        ACCEPTREJECT_PISTAR(t) = Rcpp::as<double>(output[3]);
        PIs.col(t)     = pis;
        //ALPHAS.col(t)  = alpha*vv;
    }
    return Rcpp::List::create(Rcpp::Named("PISTAR",PISTAR),
                              Rcpp::Named("RSTAR",RSTAR),
                              Rcpp::Named("PIs",PIs),
                              Rcpp::Named("ACCEPTREJECT_PISTAR",ACCEPTREJECT_PISTAR),
                              Rcpp::Named("ACCEPTREJECT_RSTAR",ACCEPTREJECT_RSTAR)
                              //Rcpp::Named("ALPHAS",ALPHAS)
    );
}
