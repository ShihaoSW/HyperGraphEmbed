# beta = 0

coverage_result = data.frame(coverage = rep(0,(400+800+1500+3000+6000)* 5 + 400 + 
                                              800 + 1500 + 3000+ 6000), 
                             length = rep(0,(400+800+1500+3000+6000)* 5 + 400 + 
                                            800 + 1500 + 3000+ 6000 ),
                             n = c(rep("n = 400", 400 * 5), rep("n = 800", 800 * 5),rep("n = 1500", 1500 * 5),
                                   rep("n = 3000", 3000 * 5),rep("n = 6000", 6000 * 5),
                                   rep("n = 400", 400), rep("n = 800", 800),rep("n = 1500", 1500),
                                   rep("n = 3000", 3000),rep("n = 6000", 6000)),
                             target = c(rep("alpha", 400), rep("F", 800), rep("Z", 800),
                                        rep("alpha", 800), rep("F", 1600), rep("Z", 1600),
                                        rep("alpha", 1500), rep("F", 3000), rep("Z", 3000),
                                        rep("alpha", 3000), rep("F", 6000), rep("Z", 6000),
                                        rep("alpha", 6000), rep("F", 12000), rep("Z", 12000),
                                        rep("Theta", 400 + 800 + 1500 + 3000+ 6000)))
coverage_result$target = factor(coverage_result$target, levels =c("alpha", "Z", "F", "Theta") )
coverage_result$n = factor(coverage_result$n, levels =c("n = 400","n = 800","n = 1500","n = 3000","n = 6000") )

p_coverage_result = data.frame(coverage = rep(0, (400 + 800 + 1500 + 3000+ 6000)*2), 
                             length = rep(0, (400 + 800 + 1500 + 3000+ 6000)*2 ),
                             n = c(rep("n = 400", 400 ), rep("n = 800", 800),rep("n = 1500", 1500),
                                   rep("n = 3000", 3000 ),rep("n = 6000", 6000 ), 
                                   rep("n = 400", 400 ), rep("n = 800", 800),rep("n = 1500", 1500),
                                   rep("n = 3000", 3000 ),rep("n = 6000", 6000 )),
                             sparsity = c(rep("0", 400 + 800 + 1500 + 3000+ 6000),
                                          rep("-1", 400 + 800 + 1500 + 3000+ 6000)))
coverage_result$target = factor(coverage_result$target, levels =c("alpha", "Z", "F", "Theta") )
coverage_result$n = factor(coverage_result$n, levels =c("n = 400","n = 800","n = 1500","n = 3000","n = 6000") )
p_coverage_result$n = factor(p_coverage_result$n, levels =c("n = 400","n = 800","n = 1500","n = 3000","n = 6000") )
p_coverage_result$sparsity = factor(p_coverage_result$sparsity, levels = c('0','-1') )

head(coverage_result)

# library(ggplot2)
# p <- ggplot(coverage_result, aes(x=target, y=coverage)) + 
#   geom_boxplot()
# p

target_cov = 0.95

n_list = c(400,800,1500,3000,6000)



for (repp in 1:1000) {
  cat(repp,"\n")
  load(paste0("./result0/result", repp, ".RData")) # change this path to where you save the results
  mean_vector = rep(0,58500 + 400 + 800 + 1500 + 3000+ 6000)
  mean_vector[1:400] = est_params[[1]]$a_mean
  mean_vector[401:800] = est_params[[1]]$F_mean[,1]
  mean_vector[801:1200] = est_params[[1]]$F_mean[,2]
  mean_vector[1201:1600] = est_params[[1]]$Z_mean[,1]
  mean_vector[1601:2000] = est_params[[1]]$Z_mean[,2]
  mean_vector[2001:2800] = est_params[[2]]$a_mean
  mean_vector[2801:3600] = est_params[[2]]$F_mean[,1]
  mean_vector[3601:4400] = est_params[[2]]$F_mean[,2]
  mean_vector[4401:5200] = est_params[[2]]$Z_mean[,1]
  mean_vector[5201:6000] = est_params[[2]]$Z_mean[,2]
  mean_vector[6001:7500] = est_params[[3]]$a_mean
  mean_vector[7501:9000] = est_params[[3]]$F_mean[,1]
  mean_vector[9001:10500] = est_params[[3]]$F_mean[,2]
  mean_vector[10501:12000] = est_params[[3]]$Z_mean[,1]
  mean_vector[12001:13500] = est_params[[3]]$Z_mean[,2]
  mean_vector[13501:16500] = est_params[[4]]$a_mean
  mean_vector[16501:19500] = est_params[[4]]$F_mean[,1]
  mean_vector[19501:22500] = est_params[[4]]$F_mean[,2]
  mean_vector[22501:25500] = est_params[[4]]$Z_mean[,1]
  mean_vector[25501:28500] = est_params[[4]]$Z_mean[,2]
  mean_vector[28501:34500] = est_params[[5]]$a_mean
  mean_vector[34501:40500] = est_params[[5]]$F_mean[,1]
  mean_vector[40501:46500] = est_params[[5]]$F_mean[,2]
  mean_vector[46501:52500] = est_params[[5]]$Z_mean[,1]
  mean_vector[52501:58500] = est_params[[5]]$Z_mean[,2]
  
  mean_vector[58501:58900] = diag(est_params[[1]]$F_mean %*% t(est_params[[1]]$Z_mean) +  
                                  rep(1,400) %*% t(est_params[[1]]$a_mean) )
  mean_vector[58901:59700] = diag(est_params[[2]]$F_mean %*% t(est_params[[2]]$Z_mean) +  
                                  rep(1,800) %*% t(est_params[[2]]$a_mean) )
  mean_vector[59701:61200] = diag(est_params[[3]]$F_mean %*% t(est_params[[3]]$Z_mean) +  
                                  rep(1,1500) %*% t(est_params[[3]]$a_mean) )
  mean_vector[61201:64200] = diag(est_params[[4]]$F_mean %*% t(est_params[[4]]$Z_mean) +  
                                  rep(1,3000) %*% t(est_params[[4]]$a_mean) )
  mean_vector[64201:70200] = diag(est_params[[5]]$F_mean %*% t(est_params[[5]]$Z_mean) +  
                                  rep(1,6000) %*% t(est_params[[5]]$a_mean) )
  
  
  width = rep(0, 58500 + 400 + 800 + 1500 + 3000+ 6000)
  for (j in 1:400) {
    width[j] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 400] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 800] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 1200] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 1600] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:800) {
    width[j + 2000] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 2800] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 3600] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 4400] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 5200] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:1500) {
    width[j + 6000] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 7500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 9000] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 10500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 12000] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:3000) {
    width[j + 13500] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 16500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 19500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 22500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 25500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:6000) {
    width[j + 28500] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 34500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 40500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 46500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 52500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][2,2])
  }
  
  cat("n = 400 \n")
  
  for (i in 1:400) {
    width[58500 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[1]]$F_mean[i,],1)) %*% est_params[[1]]$var_ests$cov_za_est[[i]] %*% c(est_params[[1]]$F_mean[i,],1) +      
           t(est_params[[1]]$Z_mean[i,]) %*% est_params[[1]]$var_ests$cov_f_est[[i]] %*%  est_params[[1]]$Z_mean[i,] )
  }
  
  cat("n = 800 \n")
  
  for (i in 1:800) {
    width[58900 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[2]]$F_mean[i,],1)) %*% est_params[[2]]$var_ests$cov_za_est[[i]] %*% c(est_params[[2]]$F_mean[i,],1) +      
             t(est_params[[2]]$Z_mean[i,]) %*% est_params[[2]]$var_ests$cov_f_est[[i]] %*%  est_params[[2]]$Z_mean[i,] )
  }
  
  cat("n = 1500 \n")
  
  for (i in 1:1500) {
    width[59700 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[3]]$F_mean[i,],1)) %*% est_params[[3]]$var_ests$cov_za_est[[i]] %*% c(est_params[[3]]$F_mean[i,],1) +      
             t(est_params[[3]]$Z_mean[i,]) %*% est_params[[3]]$var_ests$cov_f_est[[i]] %*%  est_params[[3]]$Z_mean[i,] )
  }
  
  cat("n = 3000 \n")
  
  for (i in 1:3000) {
    width[61200 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[4]]$F_mean[i,],1)) %*% est_params[[4]]$var_ests$cov_za_est[[i]] %*% c(est_params[[4]]$F_mean[i,],1) +      
             t(est_params[[4]]$Z_mean[i,]) %*% est_params[[4]]$var_ests$cov_f_est[[i]] %*%  est_params[[4]]$Z_mean[i,] )
  }
  
  cat("n = 6000 \n")
  
  for (i in 1:6000) {
    width[64200 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[5]]$F_mean[i,],1)) %*% est_params[[5]]$var_ests$cov_za_est[[i]] %*% c(est_params[[5]]$F_mean[i,],1) +      
             t(est_params[[5]]$Z_mean[i,]) %*% est_params[[5]]$var_ests$cov_f_est[[i]] %*%  est_params[[5]]$Z_mean[i,] )
  }
  

  
  cat("Theta finishes \n")
  
  low = mean_vector - width
  up = mean_vector + width
  
  
  true_mean = rep(0,58500)
  true_mean[1:400] = true_params[[1]]$alpha
  true_mean[401:800] = true_params[[1]]$F[,1]
  true_mean[801:1200] = true_params[[1]]$F[,2]
  true_mean[1201:1600] = true_params[[1]]$Z[,1]
  true_mean[1601:2000] = true_params[[1]]$Z[,2]
  true_mean[2001:2800] = true_params[[2]]$alpha
  true_mean[2801:3600] = true_params[[2]]$F[,1]
  true_mean[3601:4400] = true_params[[2]]$F[,2]
  true_mean[4401:5200] = true_params[[2]]$Z[,1]
  true_mean[5201:6000] = true_params[[2]]$Z[,2]
  true_mean[6001:7500] = true_params[[3]]$alpha
  true_mean[7501:9000] = true_params[[3]]$F[,1]
  true_mean[9001:10500] = true_params[[3]]$F[,2]
  true_mean[10501:12000] = true_params[[3]]$Z[,1]
  true_mean[12001:13500] = true_params[[3]]$Z[,2]
  true_mean[13501:16500] = true_params[[4]]$alpha
  true_mean[16501:19500] = true_params[[4]]$F[,1]
  true_mean[19501:22500] = true_params[[4]]$F[,2]
  true_mean[22501:25500] = true_params[[4]]$Z[,1]
  true_mean[25501:28500] = true_params[[4]]$Z[,2]
  true_mean[28501:34500] = true_params[[5]]$alpha
  true_mean[34501:40500] = true_params[[5]]$F[,1]
  true_mean[40501:46500] = true_params[[5]]$F[,2]
  true_mean[46501:52500] = true_params[[5]]$Z[,1]
  true_mean[52501:58500] = true_params[[5]]$Z[,2]
  true_mean[58501:58900] = diag(true_params[[1]]$F %*% t(true_params[[1]]$Z) +  
                                 rep(1,400) %*% t(true_params[[1]]$alpha) )
  true_mean[58901:59700] = diag(true_params[[2]]$F %*% t(true_params[[2]]$Z) +  
                                 rep(1,800) %*% t(true_params[[2]]$alpha) )
  true_mean[59701:61200] = diag(true_params[[3]]$F %*% t(true_params[[3]]$Z) +  
                                 rep(1,1500) %*% t(true_params[[3]]$alpha) )
  true_mean[61201:64200] = diag(true_params[[4]]$F %*% t(true_params[[4]]$Z) +  
                                 rep(1,3000) %*% t(true_params[[4]]$alpha) )
  true_mean[64201:70200] = diag(true_params[[5]]$F %*% t(true_params[[5]]$Z) +  
                                 rep(1,6000) %*% t(true_params[[5]]$alpha) )
  
  
  coverage_mid = ((true_mean <= up) & (true_mean >= low)) |(( -true_mean <= up) & (-true_mean >= low))
  coverage_result$coverage = coverage_result$coverage + coverage_mid
  coverage_result$length = coverage_result$length + width
  
  cat("doing p \n")
  p_mean_vector =  invlogit(mean_vector[58501:70200])
  p_width = b_hes(mean_vector[58501:70200]) * width[58501:70200]
  p_true_mean = invlogit(true_mean[58501:70200])
  p_low = p_mean_vector - p_width
  p_up = p_mean_vector + p_width
  p_coverage_mid = ((p_true_mean <= p_up) & (p_true_mean >= p_low)) |(( -p_true_mean <= p_up) & (-p_true_mean >= p_low))
  p_coverage_result$coverage[1:11700] = p_coverage_result$coverage[1:11700] + p_coverage_mid
  p_coverage_result$length[1:11700] = p_coverage_result$length[1:11700] + p_width
  
}

coverage_result$coverage = coverage_result$coverage / 1000
coverage_result$length = coverage_result$length / 1000 * 2

p <- ggplot(coverage_result, aes(x=target, y=coverage, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab("Coverage target") + ylab("Coverage") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Empirical coverage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Num. of vertices"))
p

p <- ggplot(coverage_result, aes(x=target, y=length, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab("Coverage target") + ylab("Interval length") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Interval length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Num. of vertices"))
p

# p <- ggplot(p_coverage_result, aes(x=sparsity, y=coverage, fill = n)) + 
#   geom_boxplot() + theme_classic() + xlab("") + ylab("Interval length") +
#   theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Interval length") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   guides(fill = guide_legend(title = "Num. of vertices"))
# p

###### beta = -1 ####

coverage_result = data.frame(coverage = rep(0,(400+800+1500+3000+6000)* 5 + 400 + 
                                              800 + 1500 + 3000+ 6000), 
                             length = rep(0,(400+800+1500+3000+6000)* 5 + 400 + 
                                            800 + 1500 + 3000+ 6000 ),
                             n = c(rep("n = 400", 400 * 5), rep("n = 800", 800 * 5),rep("n = 1500", 1500 * 5),
                                   rep("n = 3000", 3000 * 5),rep("n = 6000", 6000 * 5),
                                   rep("n = 400", 400), rep("n = 800", 800),rep("n = 1500", 1500),
                                   rep("n = 3000", 3000),rep("n = 6000", 6000)),
                             target = c(rep("alpha", 400), rep("F", 800), rep("Z", 800),
                                        rep("alpha", 800), rep("F", 1600), rep("Z", 1600),
                                        rep("alpha", 1500), rep("F", 3000), rep("Z", 3000),
                                        rep("alpha", 3000), rep("F", 6000), rep("Z", 6000),
                                        rep("alpha", 6000), rep("F", 12000), rep("Z", 12000),
                                        rep("Theta", 400 + 800 + 1500 + 3000+ 6000)))
coverage_result$target = factor(coverage_result$target, levels =c("alpha", "Z", "F", "Theta") )
coverage_result$n = factor(coverage_result$n, levels =c("n = 400","n = 800","n = 1500","n = 3000","n = 6000") )

coverage_result$target = factor(coverage_result$target, levels =c("alpha", "Z", "F", "Theta") )
coverage_result$n = factor(coverage_result$n, levels =c("n = 400","n = 800","n = 1500","n = 3000","n = 6000") )


# head(coverage_result)
# 
# p <- ggplot(coverage_result, aes(x=target, y=coverage)) + 
#   geom_boxplot()
# p

target_cov = 0.95

n_list = c(400,800,1500,3000,6000)


# setdiff(1:1000, suc)


for (repp in 1:1000) {
  cat(repp,"\n")
  load(paste0("./result1/result", repp, ".RData"))
  mean_vector = rep(0,58500 + 400 + 800 + 1500 + 3000+ 6000)
  mean_vector[1:400] = est_params[[1]]$a_mean
  mean_vector[401:800] = est_params[[1]]$F_mean[,1]
  mean_vector[801:1200] = est_params[[1]]$F_mean[,2]
  mean_vector[1201:1600] = est_params[[1]]$Z_mean[,1]
  mean_vector[1601:2000] = est_params[[1]]$Z_mean[,2]
  mean_vector[2001:2800] = est_params[[2]]$a_mean
  mean_vector[2801:3600] = est_params[[2]]$F_mean[,1]
  mean_vector[3601:4400] = est_params[[2]]$F_mean[,2]
  mean_vector[4401:5200] = est_params[[2]]$Z_mean[,1]
  mean_vector[5201:6000] = est_params[[2]]$Z_mean[,2]
  mean_vector[6001:7500] = est_params[[3]]$a_mean
  mean_vector[7501:9000] = est_params[[3]]$F_mean[,1]
  mean_vector[9001:10500] = est_params[[3]]$F_mean[,2]
  mean_vector[10501:12000] = est_params[[3]]$Z_mean[,1]
  mean_vector[12001:13500] = est_params[[3]]$Z_mean[,2]
  mean_vector[13501:16500] = est_params[[4]]$a_mean
  mean_vector[16501:19500] = est_params[[4]]$F_mean[,1]
  mean_vector[19501:22500] = est_params[[4]]$F_mean[,2]
  mean_vector[22501:25500] = est_params[[4]]$Z_mean[,1]
  mean_vector[25501:28500] = est_params[[4]]$Z_mean[,2]
  mean_vector[28501:34500] = est_params[[5]]$a_mean
  mean_vector[34501:40500] = est_params[[5]]$F_mean[,1]
  mean_vector[40501:46500] = est_params[[5]]$F_mean[,2]
  mean_vector[46501:52500] = est_params[[5]]$Z_mean[,1]
  mean_vector[52501:58500] = est_params[[5]]$Z_mean[,2]
  
  mean_vector[58501:58900] = diag(est_params[[1]]$F_mean %*% t(est_params[[1]]$Z_mean) +  
                                    rep(1,400) %*% t(est_params[[1]]$a_mean) )
  mean_vector[58901:59700] = diag(est_params[[2]]$F_mean %*% t(est_params[[2]]$Z_mean) +  
                                    rep(1,800) %*% t(est_params[[2]]$a_mean) )
  mean_vector[59701:61200] = diag(est_params[[3]]$F_mean %*% t(est_params[[3]]$Z_mean) +  
                                    rep(1,1500) %*% t(est_params[[3]]$a_mean) )
  mean_vector[61201:64200] = diag(est_params[[4]]$F_mean %*% t(est_params[[4]]$Z_mean) +  
                                    rep(1,3000) %*% t(est_params[[4]]$a_mean) )
  mean_vector[64201:70200] = diag(est_params[[5]]$F_mean %*% t(est_params[[5]]$Z_mean) +  
                                    rep(1,6000) %*% t(est_params[[5]]$a_mean) )
  
  
  width = rep(0, 58500 + 400 + 800 + 1500 + 3000+ 6000)
  for (j in 1:400) {
    width[j] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 400] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 800] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 1200] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 1600] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[1]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:800) {
    width[j + 2000] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 2800] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 3600] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 4400] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 5200] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[2]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:1500) {
    width[j + 6000] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 7500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 9000] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 10500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 12000] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[3]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:3000) {
    width[j + 13500] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 16500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 19500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 22500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 25500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[4]]$var_ests$cov_za_est[[j]][2,2])
  }
  for (j in 1:6000) {
    width[j + 28500] =   qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][3,3])
    width[j + 34500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_f_est[[j]][1,1])
    width[j + 40500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_f_est[[j]][2,2])
    width[j + 46500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][1,1])
    width[j + 52500] = qnorm( 1 - (1 - target_cov) / 2  ) * sqrt(est_params[[5]]$var_ests$cov_za_est[[j]][2,2])
  }
  
  cat("n = 400 \n")
  
  for (i in 1:400) {
    width[58500 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[1]]$F_mean[i,],1)) %*% est_params[[1]]$var_ests$cov_za_est[[i]] %*% c(est_params[[1]]$F_mean[i,],1) +      
             t(est_params[[1]]$Z_mean[i,]) %*% est_params[[1]]$var_ests$cov_f_est[[i]] %*%  est_params[[1]]$Z_mean[i,] )
  }
  
  cat("n = 800 \n")
  
  for (i in 1:800) {
    width[58900 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[2]]$F_mean[i,],1)) %*% est_params[[2]]$var_ests$cov_za_est[[i]] %*% c(est_params[[2]]$F_mean[i,],1) +      
             t(est_params[[2]]$Z_mean[i,]) %*% est_params[[2]]$var_ests$cov_f_est[[i]] %*%  est_params[[2]]$Z_mean[i,] )
  }
  
  cat("n = 1500 \n")
  
  for (i in 1:1500) {
    width[59700 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[3]]$F_mean[i,],1)) %*% est_params[[3]]$var_ests$cov_za_est[[i]] %*% c(est_params[[3]]$F_mean[i,],1) +      
             t(est_params[[3]]$Z_mean[i,]) %*% est_params[[3]]$var_ests$cov_f_est[[i]] %*%  est_params[[3]]$Z_mean[i,] )
  }
  
  cat("n = 3000 \n")
  
  for (i in 1:3000) {
    width[61200 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[4]]$F_mean[i,],1)) %*% est_params[[4]]$var_ests$cov_za_est[[i]] %*% c(est_params[[4]]$F_mean[i,],1) +      
             t(est_params[[4]]$Z_mean[i,]) %*% est_params[[4]]$var_ests$cov_f_est[[i]] %*%  est_params[[4]]$Z_mean[i,] )
  }
  
  cat("n = 6000 \n")
  
  for (i in 1:6000) {
    width[64200 + i] =  qnorm( 1 - (1 - target_cov) / 2  ) * 
      sqrt(t(c(est_params[[5]]$F_mean[i,],1)) %*% est_params[[5]]$var_ests$cov_za_est[[i]] %*% c(est_params[[5]]$F_mean[i,],1) +      
             t(est_params[[5]]$Z_mean[i,]) %*% est_params[[5]]$var_ests$cov_f_est[[i]] %*%  est_params[[5]]$Z_mean[i,] )
  }
  
  
  
  cat("Theta finishes \n")
  
  low = mean_vector - width
  up = mean_vector + width
  
  
  true_mean = rep(0,58500)
  true_mean[1:400] = true_params[[1]]$alpha
  true_mean[401:800] = true_params[[1]]$F[,1]
  true_mean[801:1200] = true_params[[1]]$F[,2]
  true_mean[1201:1600] = true_params[[1]]$Z[,1]
  true_mean[1601:2000] = true_params[[1]]$Z[,2]
  true_mean[2001:2800] = true_params[[2]]$alpha
  true_mean[2801:3600] = true_params[[2]]$F[,1]
  true_mean[3601:4400] = true_params[[2]]$F[,2]
  true_mean[4401:5200] = true_params[[2]]$Z[,1]
  true_mean[5201:6000] = true_params[[2]]$Z[,2]
  true_mean[6001:7500] = true_params[[3]]$alpha
  true_mean[7501:9000] = true_params[[3]]$F[,1]
  true_mean[9001:10500] = true_params[[3]]$F[,2]
  true_mean[10501:12000] = true_params[[3]]$Z[,1]
  true_mean[12001:13500] = true_params[[3]]$Z[,2]
  true_mean[13501:16500] = true_params[[4]]$alpha
  true_mean[16501:19500] = true_params[[4]]$F[,1]
  true_mean[19501:22500] = true_params[[4]]$F[,2]
  true_mean[22501:25500] = true_params[[4]]$Z[,1]
  true_mean[25501:28500] = true_params[[4]]$Z[,2]
  true_mean[28501:34500] = true_params[[5]]$alpha
  true_mean[34501:40500] = true_params[[5]]$F[,1]
  true_mean[40501:46500] = true_params[[5]]$F[,2]
  true_mean[46501:52500] = true_params[[5]]$Z[,1]
  true_mean[52501:58500] = true_params[[5]]$Z[,2]
  true_mean[58501:58900] = diag(true_params[[1]]$F %*% t(true_params[[1]]$Z) +  
                                  rep(1,400) %*% t(true_params[[1]]$alpha) )
  true_mean[58901:59700] = diag(true_params[[2]]$F %*% t(true_params[[2]]$Z) +  
                                  rep(1,800) %*% t(true_params[[2]]$alpha) )
  true_mean[59701:61200] = diag(true_params[[3]]$F %*% t(true_params[[3]]$Z) +  
                                  rep(1,1500) %*% t(true_params[[3]]$alpha) )
  true_mean[61201:64200] = diag(true_params[[4]]$F %*% t(true_params[[4]]$Z) +  
                                  rep(1,3000) %*% t(true_params[[4]]$alpha) )
  true_mean[64201:70200] = diag(true_params[[5]]$F %*% t(true_params[[5]]$Z) +  
                                  rep(1,6000) %*% t(true_params[[5]]$alpha) )
  
  
  coverage_mid = ((true_mean <= up) & (true_mean >= low)) |(( -true_mean <= up) & (-true_mean >= low))
  coverage_result$coverage = coverage_result$coverage + coverage_mid
  coverage_result$length = coverage_result$length + width
  
  cat("doing p \n")
  p_mean_vector =  invlogit(mean_vector[58501:70200])
  p_width = b_hes(mean_vector[58501:70200]) * width[58501:70200]
  p_true_mean = invlogit(true_mean[58501:70200])
  p_low = p_mean_vector - p_width
  p_up = p_mean_vector + p_width
  p_coverage_mid = ((p_true_mean <= p_up) & (p_true_mean >= p_low)) |(( -p_true_mean <= p_up) & (-p_true_mean >= p_low))
  p_coverage_result$coverage[11701:23400] = p_coverage_result$coverage[11701:23400] + p_coverage_mid
  p_coverage_result$length[11701:23400] = p_coverage_result$length[11701:23400] + p_width
  
}

coverage_result$coverage = coverage_result$coverage / 1000
coverage_result$length = coverage_result$length / 1000 * 2

# coverage_result$coverage[11701:23400] = coverage_result$coverage[11701:23400] * 1000
# coverage_result$length[11701:23400] = coverage_result$length[11701:23400] * 1000 * 2

p <- ggplot(coverage_result, aes(x=target, y=coverage, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab("Coverage target") + ylab("Coverage") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Empirical coverage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Num. of vertices"))
p

p <- ggplot(coverage_result, aes(x=target, y=length, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab("Coverage target") + ylab("Interval length") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Interval length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Num. of vertices"))
p

p_coverage_result$coverage = p_coverage_result$coverage / 1000
p_coverage_result$length = p_coverage_result$length / 1000 * 2


library(latex2exp)
library(ggplot2)

p <- ggplot(p_coverage_result, aes(x=sparsity, y=coverage, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab(TeX("$\\beta_{m,n}^{*}$")) + ylab("Interval length") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Empirical coverage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Empirical coverage"))
p

p <- ggplot(p_coverage_result, aes(x=sparsity, y=length, fill = n)) + 
  geom_boxplot() + theme_classic() + xlab(TeX("$\\beta_{m,n}^{*}$")) + ylab("Interval length") +
  theme(legend.position="none") +scale_fill_brewer(palette="Blues")  + ggtitle("Interval length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Empirical coverage"))
p



# p <- ggplot(coverage_result, aes(x=target, y=length, fill = n)) + 
#   geom_boxplot() + theme_classic() + xlab("Coverage target") + ylab("Interval length") +
#   theme(legend.position="top") +scale_fill_brewer(palette="Blues")  + ggtitle("Interval length") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   guides(fill = guide_legend(title = "Num. of vertices")) + 
#   theme(legend.title=element_text(size=30), 
#         legend.text=element_text(size=30),
#         legend.key.size = unit(3, 'cm'))
# p
