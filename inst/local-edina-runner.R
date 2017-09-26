library("ecdm")

library("tmsadata")

# Load data
data("trial_matrix", package="tmsadata")

# Coerce to matrix
trial_matrix = as.matrix(trial_matrix)

run_edina = function(data, k = c(3:10), burnin = 20000, chain_length = 10000) {
    # Burn in
    burnin = 20000

    # Length of chain
    chain_length = burnin + 10000

    # Coerce to matrix
    trial_matrix = as.matrix(trial_matrix)

    num_k = length(k)

    timedata = matrix(NA, num_k, 3)

    outobj = vector('list', num_k)

    for(i in seq_along(k)) {
        k_idx = k[i]
        message("Working on k = ", k_idx)
        # Launch job
        timedata[i,] = system.time({
            outobj[[i]] = dina_Gibbs_Q(trial_matrix, K=k_idx, burnin, chain_length)
        })[1:3]

        message("Time: ",  timedata[i])

        edina_obj = outobj[[i]]
        time_info = timedata[i,]
        save(edina_obj, time_info, file=paste0("sim_data_",k_idx,".rda"))
    }

    list("timing" = timedata,
         "edina_obj" = outobj)
}

# Run EDINA
d = run_edina(trial_matrix, k = c(1:4))

edina_k_1_to_4 = d

save(edina_k_1_to_4, file="edina_k_1_to_4.rda")



sample_OR = OddsRatio(nrow(trial_matrix), ncol(trial_matrix), trial_matrix)

# Every entry in the cube, look at the proportion of generated odds ratio
# that exceed the observed.

obj = matrix(0, ncol(edina_obj$ORs), ncol(edina_obj$ORs))

for(i in seq_len(dim(edina_obj$ORs)[3])){
    message(i)
    obj = obj + as.matrix((edina_obj$ORs[,,i]) > sample_OR)*1
}


d = obj /  dim(edina_obj$ORs)[3]

mu_d = mean(d[upper.tri(d)] < 0.05 | d[upper.tri(d)] > 0.95)

## model fit
## Model parameters (item parameters)
## Q matrix - take a mean, variance covariance idea... -  visualization. describe the dist - how?
## Class pi

## Find the model, describe parameter ests.

## Inference about which K
## Describe Q matrix
## Pi vector

## low dimensional representation
