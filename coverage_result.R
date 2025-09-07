#################################################
##### high-dimensional hypergraph embedding #####
############### Simulation ######################

source("functions_limit.R")

###### set seed ########
## The following commented codes are for Monte Carlo repetitions
## To reproduce results in the paper, set seeds to be (1,2,...,1000)
## seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # This line sets seed using the task ID on the Michigan Great Lakes cluster
set.seed(seed)

##### All coverage results #####
n_list = c(400,800,1500,3000,6000)

####### alpha = 0 #######
# we only record the information for the first 10 hyperedges and vertices 
true_params = list()
est_params = list() # including mean and variance estimations


K = 2
rho = 0
beta = 0

for (i in 1:length(n_list) ) {
  cat(i,"\n")
  n = n_list[i]
  m = n
  
  data = gen_data(m, n, K, beta, rho)
  V = data$hypergraph
  params_collect = list()
  params_collect$"alpha" =  data$params$alpha
  params_collect$"F" =  data$params$F_mat
  params_collect$"Z" =  data$params$Z_mat
  true_params[[i]] = params_collect
  
  initials = usvd(V, K = 5, low_p = exp(-10), 0.9)
  
  latent_init = entry_recons(initials, K)
  F0 = latent_init$F_mat
  Z0 = latent_init$Z_mat
  alpha0 = latent_init$alpha
  
  result = am_pga(V, F0, Z0, alpha0, nT = 400)
  F_est = result$F_hat
  Z_est = result$Z_hat
  alpha_est = result$alpha_hat
  Theta_est = rep(1,m) %*% t(alpha_est) + F_est %*% t(Z_est)
  
  reconstructed = entry_recons(Theta = Theta_est, K = K)
  F_new = reconstructed$F_mat
  Z_new = reconstructed$Z_mat
  alpha_new = reconstructed$alpha
  
  est_results = list()
  est_results$"F_mean" = F_new
  est_results$"Z_mean" = Z_new
  est_results$"a_mean" = alpha_new
  
  var_ests = var_estimation(F_new, Z_new, alpha_new, index_set = 1 : n_list[i])
  est_results$"var_ests" = var_ests
  
  est_params[[i]] = est_results
}

save(true_params, est_params,
     file = paste0("./result0/result", seed, ".RData") )




####### alpha = -1 #######
# we only record the information for the first 10 hyperedges and vertices 
true_params = list()
est_params = list() # including mean and variance estimations

set.seed(seed)

K = 2
rho = 0
beta = 1

for (i in 1:length(n_list) ) {
  cat(i,"\n")
  n = n_list[i]
  m = n
  
  data = gen_data(m, n, K, beta, rho)
  V = data$hypergraph
  params_collect = list()
  params_collect$"alpha" =  data$params$alpha
  params_collect$"F" =  data$params$F_mat
  params_collect$"Z" =  data$params$Z_mat
  true_params[[i]] = params_collect
  
  initials = usvd(V, K = 5, low_p = exp(-10), 0.9)
  
  latent_init = entry_recons(initials, K)
  F0 = latent_init$F_mat
  Z0 = latent_init$Z_mat
  alpha0 = latent_init$alpha
  
  result = am_pga(V, F0, Z0, alpha0, nT = 400)
  F_est = result$F_hat
  Z_est = result$Z_hat
  alpha_est = result$alpha_hat
  Theta_est = rep(1,m) %*% t(alpha_est) + F_est %*% t(Z_est)
  
  reconstructed = entry_recons(Theta = Theta_est, K = K)
  F_new = reconstructed$F_mat
  Z_new = reconstructed$Z_mat
  alpha_new = reconstructed$alpha
  
  est_results = list()
  est_results$"F_mean" = F_new
  est_results$"Z_mean" = Z_new
  est_results$"a_mean" = alpha_new
  
  var_ests = var_estimation(F_new, Z_new, alpha_new, index_set = 1:n_list[i])
  est_results$"var_ests" = var_ests
  
  est_params[[i]] = est_results
}

save(true_params, est_params,
     file = paste0("./result1/result", seed, ".RData") )


