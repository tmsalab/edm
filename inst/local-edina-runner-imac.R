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
d = run_edina(trial_matrix, k = c(5:8))
