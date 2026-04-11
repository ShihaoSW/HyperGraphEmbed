#################################################
##### high-dimensional hypergraph embedding #####
######## Supplementary material figures #########
#################################################

# This script reproduces Figures 1--6 in the supplementary material.
# Results under different parameter settings can be obtained by simply replacing
# the values in the corresponding marked locations.

################ Figures 1-3 ####################
m = 10000
n = 1000 # replace with 2000, 3000 for Figures 2 and 3
K = 4
beta = 0.5 # vary in {-0.5, -1, ..., -5.5, -6} for subfigures
rho = 0
hypergraph = gen_data(m,n,K,beta,rho)$hypergraph
orders = hypergraph %*% rep(1,n)
library(latex2exp)
hist(
  orders,
  breaks = 100,
  main = TeX("$\\beta_{m,n}^{*} = -0.5$"),
  xlab = "Orders of hyperlinks",
  ylab = "Frequency"
)

################## Figures 4 ####################
m = 2500 # vary (m,n) = (10000,2000) and K=4 for another setup
n = 500
beta = 4
rho = 0
K = 3
lkl_true = 0
lkl_usvt = rep(0,100)
lkl_or = rep(0,100)
lkl_random = rep(0,100)
for (i in 1:100) {
  cat(i,"\n")
  data = gen_data(m, n, K, beta, rho)
  V = data$hypergraph
  params = data$params
  alpha_true = params$alpha
  F_true = params$F_mat
  Z_true = params$Z_mat
  
  lkl_true = lkl_true + Log_likelihood_ori(V, alpha_true, F_true, Z_true)
  
  # usvt
  initials = usvd(V, K = 5, low_p = exp(-10), 0.9)
  
  latent_init = entry_recons(initials, K)
  F0 = latent_init$F_mat
  Z0 = latent_init$Z_mat
  alpha0 = latent_init$alpha
  
  result = am_pga(V, F0, Z0, alpha0, nT = 100)
  lkl_usvt = lkl_usvt + result$lkl 
  
  # or
  covz = diag(K)
  for (i in 1:K) {
    for (j in 1:K) {
      covz[i, j] = rho^(abs(i - j))
    }
  }
  covz = 0.2 * covz
  
  covf = diag(K)
  for (i in 1:K) {
    for (j in 1:K) {
      covf[i, j] = 0.5^(abs(i - j))
    }
  }
  covf = 0.2 * covf
  
  covz_root = Sigma_root(covz)
  covf_root = Sigma_root(covf)
  
  alpha0 = runif(n, min = -1, max = 1)
  alpha0 = alpha0 - mean(alpha0) + beta * rep(1,n)
  
  zc = diag(K)
  n_center = K
  
  F0 = matrix(rnorm(m * K), ncol = K) %*% covf_root
  F0 = trunc_mat(F0, low = -1, high = 1)
  F0 = F0 - rep(1, m) %*% t(colMeans(F0))
  
  Z0 = matrix(0, nrow = n, ncol = K)
  size = floor(n / n_center)
  for (cent in 1:(n_center - 1)) {
    matrix_add = trunc_mat(
      matrix(rnorm(size * K), ncol = K) %*% covz_root,
      low = -1, high = 1
    )
    Z0[((cent - 1) * size + 1):(cent * size), ] =
      rep(1, size) %*% t(zc[cent, ]) + matrix_add
  }
  Z0[((n_center - 1) * size + 1):n, ] =
    rep(1, n - (n_center - 1) * size) %*% t(zc[n_center, ]) +
    matrix(rnorm((n - (n_center - 1) * size) * K), ncol = K) %*% covz_root
  
  result = am_pga(V, F0, Z0, alpha0, nT = 100)
  lkl_or = lkl_or + result$lkl 
  
  # random
  F0 = matrix(rnorm(m * K, mean = 0, sd = 1), nrow = m, ncol = K)
  Z0 = matrix(rnorm(n * K, mean = 0, sd = 1), nrow = n, ncol = K)
  alpha0 = runif(n, min = -1, max = 1)
  
  result = am_pga(V, F0, Z0, alpha0, nT = 100)
  lkl_random = lkl_random + result$lkl 
}

lkl_true = lkl_true / (100 * n)
lkl_usvt = lkl_usvt / (100 * n)
lkl_or = lkl_or / (100 * n)
lkl_random = lkl_random / (100 * n)

# plot
library(ggplot2)

df <- data.frame(
  Iterations = rep(1:100, 3),
  Likelihood = c(lkl_usvt, lkl_or, lkl_random),
  method = factor(
    rep(c("Initial–USVT", "Oracle–Random", "Random"), each = 100),
    levels = c("Initial–USVT", "Oracle–Random", "Random")
  )
)

ggplot(df, aes(x = Iterations, y = Likelihood, color = method, fill = method)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = lkl_true, linetype = "dashed", color = "purple", linewidth = 1) +
  scale_color_manual(values = c("Initial–USVT" = "tomato", "Oracle–Random" = "green3", "Random" = "cornflowerblue")) +
  scale_fill_manual(values = c("Initial–USVT" = "tomato", "Oracle–Random" = "green3", "Random" = "cornflowerblue")) +
  labs(x = "Iterations", y = "Likelihood", color = "method", fill = "method") +
  theme_gray(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_line(linewidth = 0.3)
  )


################## Figure 5 ####################
m = 200
n = 200
K = 1
beta = 0.5
set.seed(2023)

# parameters
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
covz_root = sqrt(covz)
covf_root = sqrt(covf)

alpha_vals = runif(n, min = -1, max = 1)
alpha_vals = alpha_vals - mean(alpha_vals)

zc = diag(K)

F_mat = matrix(rnorm(m*K), ncol = K) %*% covf_root 
F_mat = trunc_mat(F_mat, low = -1, high = 1)
F_mat = F_mat - rep(1,m) %*% t(colMeans(F_mat))


Z_mat = matrix(rep(zc, n*K), ncol = K) + matrix(rnorm(n *K), ncol = K) %*% covz_root

Theta = F_mat %*% t(Z_mat) + rep(1,m) %*% t(alpha_vals) - beta * rep(1,m) %*% t(rep(1,n))

g = sqrt(sqrt(t(F_mat)%*%F_mat/m)) /sqrt(sqrt(t(Z_mat)%*%Z_mat /n))


alpha_true = alpha_vals - beta
F_true = F_mat %*% solve(g)
Z_true = Z_mat %*% g
t(F_true)%*%F_true 
t(Z_true)%*%Z_true

# theoretical variance of alpha estimates, replace alpha to visualize other parameters
var_ests = var_estimation(F_true, Z_true, alpha_true, index_set = 1:10)
alpha_var = var_ests$cov_za_est[[1]][2, 2]

# emprical alpha estimates distribution
alpha_list = rep(0,2000) 
for (rep in 1:2000) {
  cat(rep, "\n")
  P = invlogit(Theta)
  V = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      V[i,j] = rbinom(1,1, P[i,j])
    }
  }
  initials = usvd(V, K = 5, low_p = exp(-10), 0.9)
  
  alpha0 = colMeans(initials)
  initials = initials - rep(1,m) %*% t(alpha0)
  initials.svd = svd(initials)
  fac = matrix(initials.svd$u[,1], ncol = 1) %*% initials.svd$d[1]
  load = matrix(initials.svd$v[,1], ncol = 1)
  
  eigs = eigen(t(fac) %*% fac / m)
  
  g = sqrt(sqrt(t(fac)%*%fac/m)) /sqrt(sqrt(t(load)%*%load /n))
  F0 = fac %*% solve(g)
  Z0 = load %*% g
  
  result = am_pga(V, F0, Z0, alpha0, nT = 20)
  alpha_hat = result$alpha_hat
  alpha_list[rep] = (alpha_hat[1] - alpha_true[1])/sqrt(alpha_var)
}


hist(
  alpha_list,
  prob = TRUE,
  breaks = 25,
  col = "white",
  border = "lightgray",
  density = 25,
  angle = 45,
  xlim = c(-4, 4),
  ylim = c(0, 0.42),
  main = TeX("Density of standardized $\\hat{\\alpha}_1$"),
  xlab = TeX("Standardized $\\hat{\\alpha}_1$"),
  ylab = TeX("Density")
)

curve(dnorm(x), from = -4, to = 4, add = TRUE, col = "blue", lwd = 2)

################## Figure 6 ####################
m = 100
n = 100
K = 1
beta = 0.5
set.seed(2023)

# parameters
covz = 1
covf = 1
covz_root = sqrt(covz)
covf_root = sqrt(covf)

alpha_vals = runif(n, min = -1, max = 1)
alpha_vals = alpha_vals - mean(alpha_vals)

zc = diag(K)

F_mat = matrix(rnorm(m*K), ncol = K) %*% covf_root 
F_mat = trunc_mat(F_mat, low = -1, high = 1)
F_mat = F_mat - rep(1,m) %*% t(colMeans(F_mat))


Z_mat = matrix(rep(zc, n*K), ncol = K) + matrix(rnorm(n *K), ncol = K) %*% covz_root

Theta = F_mat %*% t(Z_mat) + rep(1,m) %*% t(alpha_vals) - beta * rep(1,m) %*% t(rep(1,n))

g = sqrt(sqrt(t(F_mat)%*%F_mat/m)) /sqrt(sqrt(t(Z_mat)%*%Z_mat /n))


alpha_true = alpha_vals - beta
F_true = F_mat %*% solve(g)
Z_true = Z_mat %*% g
t(F_true)%*%F_true 
t(Z_true)%*%Z_true

# alpha 1:20
alpha_list = alpha_true[1:20]

# emprical alpha estimates distribution
P = invlogit(Theta)
V = matrix(0, nrow = m, ncol = n)
for (i in 1:m) {
  for (j in 1:n) {
    V[i,j] = rbinom(1,1, P[i,j])
  }
}
initials = usvd(V, K = 5, low_p = exp(-10), 0.9)

alpha0 = colMeans(initials)
initials = initials - rep(1,m) %*% t(alpha0)
initials.svd = svd(initials)
fac = matrix(initials.svd$u[,1], ncol = 1) %*% initials.svd$d[1]
load = matrix(initials.svd$v[,1], ncol = 1)

eigs = eigen(t(fac) %*% fac / m)

g = sqrt(sqrt(t(fac)%*%fac/m)) /sqrt(sqrt(t(load)%*%load /n))
F0 = fac %*% solve(g)
Z0 = load %*% g

result = am_pga(V, F0, Z0, alpha0, nT = 200)
alpha_est = result$alpha_hat
F_est = result$F_hat
Z_est = result$Z_hat
var_ests = var_estimation(F_est, Z_est, alpha_est, index_set = 1:20)
alpha_up = rep(0,20)
alpha_low = rep(0,20)
for (i in 1:20) {
  alpha_up[i] = alpha_hat[i] + sqrt(var_ests$cov_za_est[[i]][2,2])
  alpha_low[i] = alpha_hat[i] - sqrt(var_ests$cov_za_est[[i]][2,2])
}

plot(
  1:20, alpha_list,
  type = "l",
  lwd = 2,
  col = "black",
  ylim = range(c(alpha_low, alpha_list, alpha_up)),
  xlab = TeX("Index $j$"),
  ylab = "",
  main = TeX("Confidence intervals for $\\{\\alpha_j^{\\dagger}\\},\\ j \\in [20]$")
)

lines(1:20, alpha_up, col = "red", lwd = 2, lty = 2)
lines(1:20, alpha_low, col = "red", lwd = 2, lty = 2)



