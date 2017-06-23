##########################################################
# Simulation runner: RRUM                                #
#                                                        #
##########################################################
# Example Call:                                          #
#                                                        #
# Rscript icc-test-rrum.R --args 3000 3 30 "results.rda" #
#                                                        #
##########################################################

library("ecdm")

message("Running simulation with ecdm v",
        as.character(utils::packageVersion("ecdm")),
        sep = "")

# Param values
args = commandArgs(trailingOnly = TRUE)

# Assign sim helpers
N        = as.numeric(args[2])
K        = as.numeric(args[3])
J        = as.numeric(args[4])
save_loc = args[6]

# Sample true attribute profiles
Z         = matrix(rnorm(N*K), N, K)
Sig       = matrix(.5, K, K)
diag(Sig) = 1
theta     = Z%*%chol(Sig)

thvals    = matrix(qnorm((1:K)/(K+1)),
                   N, K, byrow=T)

Alphas    = 1*(theta>thvals)

# Defining matrix of possible attribute profiles
As = as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1)))
Q = rbind(As[rep(c(2, 3, 5),4),],
          As[rep(c(4, 6, 7),4),],
          As[rep(8, 6),])

a = As %*% bijectionvector(K)
As = As[a+1,]

#Setting item parameters
pistar = rep(.9, J)
rstar = matrix(.6, J, K)*Q

# Simulate data under rRUM model
Y = simrRUM(N, J, K, Q, rstar, pistar, Alphas)

# Estimation Settings
chainLength = 20000
burnin = chainLength/2 + 1

# Gibbs Estimation
gibbs_estimation = rRUM_Gibbs(Y, Q, chainLength)

# Metropolis-Hastings
# mh_estimation = rRUM_MH(Y, Q, .005, chainLength)

# Release
save(gibbs_estimation,
     file = save_loc)
