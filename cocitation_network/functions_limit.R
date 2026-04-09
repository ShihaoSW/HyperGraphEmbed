###############################################
######## functions for h-d hypergraphs ########
###############################################

##### generate data and parameters #####
trunc_mat <- function(A, low, high){
  m = dim(A)[1]
  n = dim(A)[2]
  for (i in 1:m) {
    for (j in 1:n) {
      if(A[i,j]<low){A[i,j] = low}
      if(A[i,j]>high){A[i,j] = high}
    }
  }
  return(A)
}

gen_data = function(m, n, K, beta, rho, trunc = 1, alpha_bound = 1){
  
  ## generate a h-d sparse hypergraph
  ## Input: m: number of edges
  ##        n: number of nodes
  ##        K: latent space dimension
  ##        beta: the sparsity adjusted parameter
  ##        rho: correlation among latent dimensions
  ##        trunc: truncation length (half) in truncated Gaussian distributions
  ##        alpha_bound: alphas are centered i.i.d.  Uniform[-alpha_bound,alpha_bound]
  
  # covariance matrix
  covz = diag(K)
  for (i in 1:K) {
    for (j in 1:K) {
      covz[i,j] = rho^( abs(i-j) )
    }
  }
  covz = 0.2 * covz
  covf = diag(K)
  for (i in 1:K) {
    for (j in 1:K) {
      covf[i,j] = 0.5^( abs(i-j) )
    }
  }
  covf = 0.2 * covf
  covz_root = Sigma_root(covz)
  covf_root = Sigma_root(covf)
  
  alpha_vals = runif(n, min = -alpha_bound, max = alpha_bound)
  alpha_vals = alpha_vals - mean(alpha_vals)
  
  # Natural parameters
  zc = diag(K)
  n_center = K
  
  F_mat = matrix(rnorm(m*K), ncol = K) %*% covf_root 
  F_mat = trunc_mat(F_mat, low = -trunc, high = trunc)
  F_mat = F_mat - rep(1,m) %*% t(colMeans(F_mat))
  
  Z_mat = matrix(0, ncol = K, nrow = n)
  size = floor(n / n_center)
  for (cent in 1: (n_center - 1)) {
    matrix_add = trunc_mat(matrix(rnorm(size *K), ncol = K) %*% covz_root, low = -trunc, high = trunc )
    Z_mat[((cent-1)*size + 1): (cent * size) , ] = rep(1, size) %*% t(zc[cent,]) + matrix_add
  }
  Z_mat[((n_center-1)*size + 1): n , ] = rep(1, (n -  (n_center-1)*size ))  %*%  t(zc[n_center,]) +
    matrix(rnorm((n -  (n_center-1)*size ) *K), ncol = K) %*% covz_root
  
  Theta = F_mat %*% t(Z_mat) + rep(1,m) %*% t(alpha_vals) - beta * rep(1,m) %*% t(rep(1,n))
  
  params = entry_recons(Theta, K)

  # Probability matrix
  P = invlogit(Theta)
  
  # Hyper edges
  hypergraph = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      hypergraph[i,j] = rbinom(1,1, P[i,j])
    }
  }
  
  res = list(hypergraph = hypergraph, params = params) 
  
}

##### logit, inv, b'', norms, likelihood #####
invlogit <- function(x){( 1 / (1 + exp(-x)))}
logit <- function(x){log(x / (1 - x))}
b_hes <- function(x){ exp(x) / (1 + exp(x))^2 }
Fnorm <- function(x){ sqrt(sum(diag(t(x) %*% x))  ) }
l2norm <- function(x){ t(x)%*%x }
op_norm <- function(x){
  sigma = t(x) %*% x
  eig_res = eigen(sigma)
  return(eig_res$values[1])
}
Log_likelihood <- function(V, F0, Z0){
  m = dim(V)[1]
  n = dim(V)[2]
  Theta = F0 %*% t(Z0)
  summa = Theta * V - log(1 + exp(Theta) )
  summa = sum(summa)
  return(summa / m)
}

Log_likelihood_ori <- function(V, alpha, F0, Z0){
  m = dim(V)[1]
  n = dim(V)[2]
  Theta = rep(1,m) %*% t(alpha) + F0 %*% t(Z0)
  summa = Theta * V - log(1 + exp(Theta) )
  summa = sum(summa)
  return(summa / m)
}

##### trim the matrix #####
trim_Theta <- function(Theta, low, up){
  m = dim(Theta)[1]
  n = dim(Theta)[2]
  for (i in 1:m) {
    for (j in 1:n) {
      if(Theta[i,j] < low){Theta[i,j] = low}
      if(Theta[i,j] > up){Theta[i,j] = up}
    }
  }
  return(Theta)
}

trim_p <-function(P, low, up){
  m = dim(P)[1]
  n = dim(P)[2]
  for (i in 1:m) {
    for (j in 1:n) {
      if(P[i,j] < low){P[i,j] = low}
      if(P[i,j] > up){P[i,j] = up}
    }
  }
  return(P)
}

trim_vec <- function(x, low, up){
  l = length(x)
  for (i in 1:l) {
    if(x[i] < low){x[i] = low}
    if(x[i] > up){x[i] = up}
  }
  return(x)
}



##### reconstruct from Theta #####
ori_int_recons <- function(Theta, K){
  alpha = colMeans(Theta)
  m = dim(Theta)[1]
  Theta = Theta - rep(1,m) %*% t(alpha)
  Theta.svd = svd(Theta)
  fac = Theta.svd$u[,1:K] %*% diag(Theta.svd$d[1:K]) / sqrt(n)
  load = Theta.svd$v[,1:K] * sqrt(n)
  return(list(F0 = fac, Z0 = load, alpha0 = alpha))
}

latent_recons <- function(Theta, K){
  n = dim(Theta)[2]
  Theta.svd = svd(Theta)
  fac = Theta.svd$u[,1:K] %*% diag(Theta.svd$d[1:K]) / sqrt(n)
  load = Theta.svd$v[,1:K] * sqrt(n)
  return(list(Factors = fac, Loadings = load))
}
entry_recons <- function(Theta, K){
  m = dim(Theta)[1]
  alpha = colMeans(Theta)
  Theta = Theta - rep(1,m) %*% t(alpha)
  Theta.svd = svd(Theta)
  fac = Theta.svd$u[,1:K] %*% diag(Theta.svd$d[1:K])
  load = Theta.svd$v[,1:K]
  
  eigs = eigen(t(fac) %*% fac / m)
  
  G = diag(1/sqrt(Theta.svd$d[1:K] ) ) * sqrt(sqrt(m/n))
  #G = eigs$vectors %*% diag(1 / sqrt(sqrt(eigs$values)) )
  Z_new = load %*% t(solve(G))
  F_new = fac %*% G
  return(list(alpha = alpha, F_mat = F_new, Z_mat = Z_new))
}

##### generate the hypergraph data #####
Sigma_root <- function(Sigma){
  eig_res = eigen(Sigma) 
  eig_space = eig_res$vectors
  eig_val = eig_res$values
  return(eig_space %*% diag(sqrt(eig_val)) %*% t(eig_space) )
}

gen_graph_ori <- function(m, n, zc, betas, alphas = c(-1,1), covz, covf){
  
  ## generate a h-d sparse hypergraph
  ## Input: m: number of edges
  ##        n: number of nodes
  ##        zc: center of the loadings
  ##        betas: beta star, the sparsity adjusted parameter
  ##        covz: covariance of the loading vector
  ##        covf: covariance of the factor vector
  ##        theta_low: lower bound for the natural parameters
  ##        theta_high: upper bound for the natural parameters
  
  covz_root = Sigma_root(covz)
  covf_root = Sigma_root(covf)
  
  alpha_vals = runif(n, min = alphas[1], max = alphas[2])
  alpha_vals = alpha_vals - mean(alpha_vals)
  
  # Natural parameters
  K = dim(zc)[2]
  n_center = dim(zc)[1]
  
  F_mat = matrix(rnorm(m*K), ncol = K) %*% covf_root 
  F_mat = F_mat - rep(1,m) %*% t(colMeans(F_mat))
  
  Z_mat = matrix(0, ncol = K, nrow = n)
  size = floor(n / n_center)
  for (cent in 1: (n_center - 1)) {
    Z_mat[((cent-1)*size + 1): (cent * size) , ] = rep(1, size) %*% t(zc[cent,]) + matrix(rnorm(size *K), ncol = K) %*% covz_root 
  }
  Z_mat[((n_center-1)*size + 1): n , ] = rep(1, (n -  (n_center-1)*size ))  %*%  t(zc[n_center,]) +
    matrix(rnorm((n -  (n_center-1)*size ) *K), ncol = K) %*% covz_root
  
  Theta = F_mat %*% t(Z_mat) + rep(1,m) %*% t(alpha_vals) - betas * rep(1,m) %*% t(rep(1,n))
  
  # Probability matrix
  P = invlogit(Theta)
  
  # Hyper edges
  hypergraph = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      hypergraph[i,j] = rbinom(1,1, P[i,j])
    }
  }
  
  V = hypergraph
  
  # trim null edges and nodes
  trim_F = c()
  trim_Z = c()
  for (i in 1:m) {
    if(length(which(V[i,] == 1)) <= 1){
      trim_F = c(trim_F, i)
    }
  }
  for (j in 1:n) {
    if(length(which(V[,j] == 1)) == 0){
      trim_Z = c(trim_Z, j)
    }
  }
  
  rem_F = setdiff(1:m, trim_F)
  rem_Z = setdiff(1:n, trim_Z)
  
  F_mat = F_mat[rem_F,]
  Z_mat = Z_mat[rem_Z,]
  alpha_vals = alpha_vals[rem_Z]
  V = V[rem_F, ]
  V = V[,rem_Z]
  
  res = list(graph = V, F_true = F_mat, Z_true = Z_mat, alpha_true = alpha_vals - betas) 
  return(res)
}

gen_graph <- function(m, n, zc, betas, covz, covf){
  
  ## generate a h-d sparse hypergraph
  ## Input: m: number of edges
  ##        n: number of nodes
  ##        zc: center of the loadings
  ##        betas: beta star, the sparsity adjusted parameter
  ##        covz: covariance of the loading vector
  ##        covf: covariance of the factor vector
  ##        theta_low: lower bound for the natural parameters
  ##        theta_high: upper bound for the natural parameters
  
  covz_root = Sigma_root(covz)
  covf_root = Sigma_root(covf)
  
  # Natural parameters
  K = dim(zc)[2]
  n_center = dim(zc)[1]
  
  F_mat = - betas * rep(1,m) %*% t(rep(1, K)) +  matrix(rnorm(m*K), ncol = K) %*% covf_root 
  
  Z_mat = matrix(0, ncol = K, nrow = n)
  size = floor(n / n_center)
  for (cent in 1: (n_center - 1)) {
    Z_mat[((cent-1)*size + 1): (cent * size) , ] = rep(1, size) %*% t(zc[cent,]) + matrix(rnorm(size *K), ncol = K) %*% covz_root 
  }
  Z_mat[((n_center-1)*size + 1): n , ] = rep(1, (n -  (n_center-1)*size ))  %*%  t(zc[n_center,]) +
    matrix(rnorm((n -  (n_center-1)*size ) *K), ncol = K) %*% covz_root
  
  Theta = F_mat %*% t(Z_mat) 
  
  # Probability matrix
  P = invlogit(Theta)
  
  # Hyper edges
  hypergraph = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      hypergraph[i,j] = rbinom(1,1, P[i,j])
    }
  }
  
  V = hypergraph
  
  # trim null edges and nodes
  trim_F = c()
  trim_Z = c()
  for (i in 1:m) {
    if(length(which(V[i,] == 1)) <= 1){
      trim_F = c(trim_F, i)
    }
  }
  for (j in 1:n) {
    if(length(which(V[,j] == 1)) == 0){
      trim_Z = c(trim_Z, j)
    }
  }
  
  rem_F = setdiff(1:m, trim_F)
  rem_Z = setdiff(1:n, trim_Z)
  
  latent_recon = latent_recons(Theta, K)
  F_mat = latent_recon$Factors
  Z_mat = latent_recon$Loadings 
  
  F_mat = F_mat[rem_F,]
  Z_mat = Z_mat[rem_Z,]
  V = V[rem_F, ]
  V = V[,rem_Z]
  
  res = list(graph = V, F_true = F_mat, Z_true = Z_mat) 
  return(res)
}

##### Initialization via USVD #####
usvd <- function(V, K, low_p, up_p){
  
  ## Universal singular value decomposition for hyper-graphs
  ## Input: V: hyper-dges
  ##        K: rank of svd
  ##        low_p: lower bound for probability matrix, exp(-beta)
  ##        up_p: upper bound for probability matrix, 1/2
  
  V.svd = svd(V)
  sigmas = V.svd$d[1:K]
  u = V.svd$u[,1:K]
  v = V.svd$v[,1:K]
  P0 = u %*% diag(sigmas) %*% t(v)
  
  P = trim_p(P0, low_p, up_p)
  Theta = logit(P)
  
  return(Theta)
} 

##### Searching for the estimator #####
am_pga <- function(Y, F0, Z0, alpha0, eta = 1, nT = 100){
  ## Projected gradient ascent for embedding models 
  ## Input: Y: the hyperedges
  ##        F0: initial values for factors
  ##        Z0: initial values for loadings
  ##        eta: start step parameter in gradient descent
  ##        nT: number of iterations
  
  # step sizes
  m = dim(Y)[1]
  n = dim(Y)[2]
  op_F = op_norm(F0)
  op_Z = op_norm(Z0)
  eta_F = 2 * eta / (op_F + op_Z)
  eta_Z = 2 * eta / (op_F + op_Z)
  eta_alpha = eta / (2 * m)
  
  F_new = F0
  Z_new = Z0
  alpha_new = alpha0 
  
  lkl = rep(0, nT)
  lkl[1] = Log_likelihood_ori(V, alpha0, F0, Z0)
  
  # iterates
  for (t in 2:nT) {
    F_old = F_new
    Z_old = Z_new
    alpha_old = alpha_new
    Theta = rep(1,m) %*% t(alpha_old) + F_old %*% t(Z_old)
    F_new = F_old + eta_F * (Y - invlogit(Theta)) %*% Z_old
    Z_new = Z_old + eta_Z * t(Y - invlogit(Theta)) %*% F_old
    alpha_new = alpha_old + eta_alpha * t(Y - invlogit(Theta)) %*% rep(1,m)
    F_new = F_new - rep(1,m) %*% t(colMeans(F_new))
    lkl[t] = Log_likelihood_ori(Y, alpha_new, F_new, Z_new)
    cat("The",t,"th iteration: ", Log_likelihood_ori(Y, alpha_new, F_new, Z_new),"\n")
  }
  
  return(list(F_hat = F_new, Z_hat = Z_new, alpha_hat = alpha_new, lkl = lkl))
  
}

##### variance estimation #####
var_estimation <- function(F_est, Z_est, alpha_est, index_set = 1:10){
  m = dim(F_est)[1]
  n = dim(Z_est)[1]
  
  Theta_hat = F_est %*% t(Z_est) + rep(1,m) %*% t(alpha_est)
  hes_hat = b_hes(Theta_hat)
  F_aug_hat = cbind(F_est, rep(1,m))
  K = dim(F_est)[2]
  
  cov_za_est = list()
  for (j in index_set) {
    omega = matrix(0, nrow = (K+1), ncol = (K+1))
    for (i in 1:m) {
      omega = omega + hes_hat[i,j] * F_aug_hat[i,] %*% t(F_aug_hat[i,])
    }
    cov_za_est[[j]] = solve(omega)
  }
  
  cov_f_est = list()
  for (i in index_set) {
    omega = matrix(0, nrow = K, ncol = K)
    for (j in 1:n) {
      omega = omega + hes_hat[i,j] * Z_est[j,] %*% t(Z_est[j,])
    }
    cov_f_est[[i]] = solve(omega)
  }
  
  
  return(list(cov_za_est = cov_za_est, cov_f_est = cov_f_est ))  
}



