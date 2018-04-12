#' Exploratory reduced Reparameterized Unified Model (ErRUM)
#'
#' Obtains samples from posterior distributon for the Exploratory
#' reduced Reparametrized Unified Model (rRUM).
#' @inheritParams edina
#'
#' @export
#' @return An `errum` object that has:
#' - `PISTAR`
#' - `RSTAR`
#' - `PIs`
#' - `QS`
#' - `m_Delta`
#' - `Delta_biject`
#' - `M2`
#' - `M1`
#' - `NUS`
#' @examples
#' \dontrun{
#' # Assign sim helpers
#' N        = 3000
#' K        = 3
#' J        = 30
#'
#' # Sample true attribute profiles
#' Z         = matrix(rnorm(N*K), N, K)
#' Sig       = matrix(.5, K, K)
#' diag(Sig) = 1
#' theta     = Z%*%chol(Sig)
#'
#' thvals    = matrix(qnorm((1:K)/(K+1)),
#'                    N, K, byrow=T)
#'
#' Alphas    = 1*(theta > thvals)
#'
#' # Defining matrix of possible attribute profiles
#' As = as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1)))
#' Q = rbind(As[rep(c(2, 3, 5),4),],
#'           As[rep(c(4, 6, 7),4),],
#'           As[rep(8, 6),])
#'
#' a = As %*% bijectionvector(K)
#' As = As[a+1,]
#'
#' # Setting item parameters
#' pistar = rep(.9, J)
#' rstar = matrix(.6, J, K)*Q
#'
#' # Simulate data under rRUM model
#' Y = sim_rrum(Q, rstar, pistar, Alphas)
#'
#' # Estimation Settings
#' chainLength = 20000
#' burnin = chainLength/2 + 1
#'
#' # Gibbs Estimation
#' gibbs_estimation = errum(Y, K, burnin, chainLength)
#' }
errum = function(data, k = 3,
                 burnin = 10000, chain_length = 10000) {
    o = errum_Gibbs_Q(data, k, burnin, chain_length)
    class(o) = c("errum")
    o
}
