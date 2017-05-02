options(width = 10000, max.print = 100000)

# Set library via bash profile
# ~/Rlib
library("ecdm")

cat("Running simulation with ecdm v", utils::packageVersion("ecdm"), "\n")

# Parm values
cmdargs = as.numeric(commandArgs(trailingOnly = TRUE))
rep     = cmdargs[1]
N       = cmdargs[2]
rho     = cmdargs[3]
K       = cmdargs[4]

compval = cmdargs[5]

J = 18


# Specify Q
if (K == 3) {
    qbj = c(4, 2, 1, 4, 2, 1, 4, 2, 1, 6, 5, 3, 6, 5, 3, 7, 7, 7)
}else if (K == 4) {
    qbj = c(8, 4, 2, 1, 8, 4, 2, 1, 12, 10, 9, 6, 5, 3, 14, 13, 11, 7)
}

# Fill Q Matrix
Q = matrix(, J, K)
for (j in 1:J) {
    Q[j, ] = inv_bijectionvector(K, qbj[j])
}

# Item parm vals
ss = gs = rep(.2, J)

# Generating attribute classes
if (rho == 0) {
    PIs = rep(1 / (2 ^ K), 2 ^ K)
    CLs = c((1:(2 ^ K)) %*% rmultinom(n = N, size = 1, prob = PIs)) - 1
}

if (rho > 0) {
    Z = matrix(rnorm(N * K), N, K)
    Sig = matrix(rho, K, K)
    diag(Sig) = 1
    X = Z %*% chol(Sig)
    thvals = matrix(rep(0, K), N, K, byrow = T)
    Alphas = 1 * (X > thvals)
    CLs = Alphas %*% bijectionvector(K)
}

# Simulate data under DINA model
ETA = ETAmat(K, J, Q)
Y_sim = sim_Y_dina(N, J, CLs, ETA, gs, ss)

# Estimating Q
vj = bijectionvector(J)
Qvj = t(Q) %*% vj
burnin = 20000
chain_length = burnin + 10000
system.time(out <- dina_Gibbs_Q(Y_sim, K, burnin, chain_length))

# Q matrix
m_Qs = apply(out$QS, c(1, 2), mean)
m_Qsvj = t(m_Qs) %*% vj
# m_Qsvj
cbind(Q[, order(Qvj)], m_Qs[, order(m_Qsvj)])

Qest = 1 * (m_Qs > .5)
Qest_check = (Q[, order(Qvj)] == Qest[, order(m_Qsvj)])

# min(Qest_check)
# mean(Qest_check)
# apply(Qest_check, 1, min)

# Parameters
m_gs = apply(out$GS, 1, mean)
# m_gs
m_ss = apply(out$SS, 1, mean)
# m_ss
m_pi = apply(out$PIs, 1, mean)
# m_pi

# Save results
save(m_gs, m_ss, m_pi, Qest, Qest_check, m_Qsvj, file = compval)
