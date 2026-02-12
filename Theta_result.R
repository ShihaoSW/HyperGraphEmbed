#################################################
##### high-dimensional hypergraph embedding #####
############### Simulation ######################

source("functions_limit.R")

###### set seed ########
## The following commented codes are for Monte Carlo repetitions
## To reproduce results in the paper, set seeds to be (1,2,...,100)
## seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # This line sets seed using the task ID on the Michigan Great Lakes cluster
set.seed(seed)



##### Fix beta = -3, change n, rho = 0 #####
n_list = seq(100,1000,100)
K_list = c(2,3,4)

errors1 = list()

rho = 0
beta = -3

for (i in 1:length(n_list) ) {
  
  error_list = c(0,0,0)
  
  for (j in 1:length(K_list)) {
    
    cat(i,j,"\n")
    n = n_list[i]
    K = K_list[j]
    m = 10*n
    
    data = gen_data(m, n, K, beta, rho)
    V = data$hypergraph

    alpha_t =  data$params$alpha
    F_t =  data$params$F_mat
    Z_t =  data$params$Z_mat
    Theta_t = F_t %*% t(Z_t) + rep(1,m) %*% t(alpha_t)
    
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
   
    error_list[j] = sqrt(sum((Theta_est - Theta_t) * (Theta_est - Theta_t))) / (sqrt(m * n))
    
  }
  
  errors1[[i]] = error_list
}



##### Fix beta = -3, change n, rho = 0.5 #####
n_list = seq(100,1000,100)
K_list = c(2,3,4)

errors2 = list()

rho = 0.5
beta = -3

for (i in 1:length(n_list) ) {
  
  error_list = c(0,0,0)
  
  for (j in 1:length(K_list)) {
    
    cat(i,j,"\n")
    n = n_list[i]
    K = K_list[j]
    m = 10*n
    
    data = gen_data(m, n, K, beta, rho)
    V = data$hypergraph
    
    alpha_t =  data$params$alpha
    F_t =  data$params$F_mat
    Z_t =  data$params$Z_mat
    Theta_t = F_t %*% t(Z_t) + rep(1,m) %*% t(alpha_t)
    
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
    
    error_list[j] = sqrt(sum((Theta_est - Theta_t) * (Theta_est - Theta_t))) / (sqrt(m * n))
    
  }
  
  errors2[[i]] = error_list
}



##### Fix m = 10n = 5000, change beta, rho = 0 #####
beta_list = seq(-0.5,-4,-0.5)
K_list = c(2,3,4)

errors3 = list()

rho = 0
n = 500
m = 10*n

for (i in 1:length(beta_list) ) {
  
  error_list = c(0,0,0)
  
  for (j in 1:length(K_list)) {
    
    cat(i,j,"\n")
    beta = beta_list[i]
    K = K_list[j]
    
    data = gen_data(m, n, K, beta, rho)
    V = data$hypergraph
    
    alpha_t =  data$params$alpha
    F_t =  data$params$F_mat
    Z_t =  data$params$Z_mat
    Theta_t = F_t %*% t(Z_t) + rep(1,m) %*% t(alpha_t)
    
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
    
    error_list[j] = sqrt(sum((Theta_est - Theta_t) * (Theta_est - Theta_t))) / (sqrt(m * n))
    
  }
  
  errors3[[i]] = error_list
}



##### Fix m = 10n = 5000, change beta, rho = 0.5 #####
beta_list = seq(-0.5,-4,-0.5)
K_list = c(2,3,4)

errors4 = list()

rho = 0.5
n = 500
m = 10*n

for (i in 1:length(beta_list) ) {
  
  error_list = c(0,0,0)
  
  for (j in 1:length(K_list)) {
    
    cat(i,j,"\n")
    beta = beta_list[i]
    K = K_list[j]
    
    data = gen_data(m, n, K, beta, rho)
    V = data$hypergraph
    
    alpha_t =  data$params$alpha
    F_t =  data$params$F_mat
    Z_t =  data$params$Z_mat
    Theta_t = F_t %*% t(Z_t) + rep(1,m) %*% t(alpha_t)
    
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
    
    error_list[j] = sqrt(sum((Theta_est - Theta_t) * (Theta_est - Theta_t))) / (sqrt(m * n))
    
  }
  
  errors4[[i]] = error_list
}



save(errors1, errors2, errors3, errors4,  file = paste0("./result/result", seed, ".RData") ) 