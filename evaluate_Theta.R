library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)


# Run the folloing codes to reproduce Figure 2 (a)
# settings
n_vals <- seq(100, 1000, by = 100)
K_vals <- c(2, 3, 4)

# storage array: [n_index, K_index, seed]
err_array <- array(NA, dim = c(length(n_vals), length(K_vals), 100))

for (seed in 1:100) {
  
  load(paste0("./result/result", seed, ".RData"))  # loads object: errors
  
  for (i in 1:length(n_vals)) {
    err_array[i, , seed] <- errors1[[i]]
  }
}

# average over seeds
err_mean <- apply(err_array, c(1, 2), mean)

df_plot <- expand.grid(
  n = n_vals,
  K = K_vals
)

df_plot$error <- as.vector(err_mean)

df_plot$K <- factor(df_plot$K)


ggplot(df_plot, aes(x = n, y = error, color = K)) +
  geom_line(linewidth = 1.4) +
  labs(
    x = "n",
    y = TeX("$\\frac{\\|\\hat{\\Theta}-\\Theta^{*}\\|_{F}}{(mn)^{1/2}}$"),
    color = "K"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    legend.key   = element_blank(),          # remove gray boxes
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2)     # thicker legend lines
    )
  )



# Run the folloing codes to reproduce Figure 2 (b)
# settings
n_vals <- seq(100, 1000, by = 100)
K_vals <- c(2, 3, 4)

# storage array: [n_index, K_index, seed]
err_array <- array(NA, dim = c(length(n_vals), length(K_vals), 100))

for (seed in 1:100) {
  
  load(paste0("./result/result", seed, ".RData"))  # loads object: errors
  
  for (i in 1:length(n_vals)) {
    err_array[i, , seed] <- errors2[[i]]
  }
}

# average over seeds
err_mean <- apply(err_array, c(1, 2), mean)

df_plot <- expand.grid(
  n = n_vals,
  K = K_vals
)

df_plot$error <- as.vector(err_mean)

df_plot$K <- factor(df_plot$K)


ggplot(df_plot, aes(x = n, y = error, color = K)) +
  geom_line(linewidth = 1.4) +
  labs(
    x = "n",
    y = TeX("$\\frac{\\|\\hat{\\Theta}-\\Theta^{*}\\|_{F}}{(mn)^{1/2}}$"),
    color = "K"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    legend.key   = element_blank(),          # remove gray boxes
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2)     # thicker legend lines
    )
  )



# Run the folloing codes to reproduce Figure 3 (a)
# settings
beta_vals = seq(-0.5,-4,-0.5)
K_vals <- c(2, 3, 4)

# storage array: [n_index, K_index, seed]
err_array <- array(NA, dim = c(length(beta_vals), length(K_vals), 100))

for (seed in 1:100) {
  
  load(paste0("./result/result", seed, ".RData"))  # loads object: errors
  
  for (i in 1:length(beta_vals)) {
    err_array[i, , seed] <- errors3[[i]]
  }
}

# average over seeds
err_mean <- apply(err_array, c(1, 2), mean)

df_plot <- expand.grid(
  beta = - beta_vals,
  K = K_vals
)

df_plot$error <- as.vector(err_mean)

df_plot$K <- factor(df_plot$K)


ggplot(df_plot, aes(x = beta, y = error, color = K)) +
  geom_line(linewidth = 1.4) +
  labs(
    x = TeX("$-\\beta_{mn}^*$"),
    y = TeX("$\\frac{\\|\\hat{\\Theta}-\\Theta^{*}\\|_{F}}{(mn)^{1/2}}$"),
    color = "K"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    legend.key   = element_blank(),          # remove gray boxes
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2)     # thicker legend lines
    )
  )



# Run the folloing codes to reproduce Figure 3 (b)
# settings
beta_vals = seq(-0.5,-4,-0.5)
K_vals <- c(2, 3, 4)

# storage array: [n_index, K_index, seed]
err_array <- array(NA, dim = c(length(beta_vals), length(K_vals), 100))

for (seed in 1:100) {
  
  load(paste0("./result/result", seed, ".RData"))  # loads object: errors
  
  for (i in 1:length(beta_vals)) {
    err_array[i, , seed] <- errors4[[i]]
  }
}

# average over seeds
err_mean <- apply(err_array, c(1, 2), mean)

df_plot <- expand.grid(
  beta = - beta_vals,
  K = K_vals
)

df_plot$error <- as.vector(err_mean)

df_plot$K <- factor(df_plot$K)


ggplot(df_plot, aes(x = beta, y = error, color = K)) +
  geom_line(linewidth = 1.4) +
  labs(
    x = TeX("$-\\beta_{mn}^*$"),
    y = TeX("$\\frac{\\|\\hat{\\Theta}-\\Theta^{*}\\|_{F}}{(mn)^{1/2}}$"),
    color = "K"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    legend.key   = element_blank(),          # remove gray boxes
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 2)     # thicker legend lines
    )
  )
