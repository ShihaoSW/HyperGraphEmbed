######################################################
########## plot author embeddings with CRs ###########
######################################################

# The code below includes some exploratory and experimental sections (commented out).
# To generate the figures in the paper, run only the uncommented sections.

load("hg1991-2000_3000_final_embeddings.RData")
library(RSpectra)
source("functions_limit.R")

##### calculate confidence regions for the selected authors #####
F_est = final_embeddings$F_em
Z_est = final_embeddings$Z_em
F_em = final_embeddings$F_em
Z_em = final_embeddings$Z_em
alpha_est = final_embeddings$alpha
m = dim(F_est)[1]
K = dim(F_est)[2]
n = dim(Z_est)[1]

# new selected authors
selected_authors = c("David Donoho", "Iain Johnstone", "Theo Gasser", "Byeong U. Park",
                     "Peter Hall", "Jianqing Fan", "Michael Martin", "Lawrence 1 Brown",
                     "Nicholas Polson", "Luis Pericchi", "Peter Green", "Luke Tierney",
                     "Raymond Carroll", "Bradley Efron", "Robert Tibshirani", "Trevor Hastie",
                     "Tze Leung Lai", "Jane-ling Wang", "Robert Kass", "Larry Wasserman", 
                     "Ross Prentice", "Lee-jen Wei", "Danyu Y. Lin", "John Oquigley",
                     "Janet Wittes", "Adrian Smith", "Alan Gelfand", "Scott Zeger", 
                     "Kung Yee Liang", "Nan Laird", "Norman Breslow", "Donald B. Rubin",
                     "Lueping Zhao", "Stuart Lipsitz", "Michael Proschan", "Joseph Lellouch")
authors3000
selected_authors_index = rep(0,36)
for (i in 1:36) {
  cat(length(which(authors3000 == selected_authors[i])), "\n")
  selected_authors_index[i] = which(authors3000 == selected_authors[i])[1]
}
selected_authors_index
authors3000[selected_authors_index]

selected_author_label = c("David Donoho", "Iain Johnstone", "Theo Gasser", "Byeong U. Park",
                          "Peter Hall", "Jianqing Fan", "Michael Martin", "Lawrence Brown",
                          "Nicholas Polson", "Luis Pericchi", "Peter Green", "Luke Tierney",
                          "Raymond Carroll", "Bradley Efron", "Robert Tibshirani", "Trevor Hastie",
                          "Tze Leung Lai", "Jane-Ling Wang", "Robert Kass", "Larry Wasserman", 
                          "Ross Prentice", "Lee-Jen Wei", "Danyu Y. Lin", "John Oquigley",
                          "Janet Wittes", "Adrian Smith", "Alan Gelfand", "Scott Zeger", 
                          "Kung Yee Liang", "Nan Laird", "Norman Breslow", "Donald B. Rubin",
                          "Lueping Zhao", "Stuart Lipsitz", "Michael Proschan", "Joseph Lellouch")

# reconstruct
Theta_est = rep(1,m) %*% t(alpha_est) + F_est %*% t(Z_est)
Theta_est = Theta_est - rep(1,m) %*% t(alpha_est)

K0 = 2
Theta.svd = svds(Theta_est, k = K0)
sigmas = Theta.svd$d[1:K0]
u = Theta.svd$u[,1:K0]
v = Theta.svd$v[,1:K0]
n = 3000
F_final = sqrt(sqrt(m / n)) * u %*% diag(sqrt(sigmas))
Z_final = sqrt(sqrt(n / m)) * v %*% diag(sqrt(sigmas))

which(abs(t(F_final) %*% F_final / m - t(Z_final) %*% Z_final / n) <= 1e-6 )

# asymptotic covariance 
Theta_hat = F_final %*% t(Z_final) + rep(1,m) %*% t(alpha_est)
abs(Theta_hat - Theta_est - rep(1,m) %*% t(alpha_est) )[1:10,1:10]
colMeans(F_final)
F_aug_hat = cbind(rep(1,m), F_final)
index_set = selected_authors
hes_hat = b_hes(Theta_hat)

cov_za_est = list()
for (j in 1:length(index_set)) {
  omega = matrix(0, nrow = (K+1), ncol = (K+1))
  cat(j,"\n")
  for (i in 1:m) {
    omega = omega + hes_hat[i,j] * F_aug_hat[i,] %*% t(F_aug_hat[i,])
  }
  cov_za_est[[j]] = solve(omega)
}
mean_list = Z_final[selected_authors_index, 1:2]
covariance_list = list()

for (j in 1:length(selected_authors)) {
  covariance_list[[j]] = cov_za_est[[j]][2:3,2:3]
}

# a function to compute the parameters in the ellipse
# ##### test the function #####
#library(ggplot2)
#library(ggforce)

#x1 = c(1,2)
#x2 = c(1,2)
#a = c(1,2)
#b = c(2,1)
#angle = c(3*pi/4, pi/2)

#ggplot() +
#geom_ellipse(aes(x0 = x1, y0 = x2, a = a, b = b, angle =angle),
#             color = 'red4', lty = 2)

#a = matrix(rnorm(4),ncol = 2 ) 
#a_eig = eigen(t(a) %*% a )

#axises = c(4,1)
#Sigma = a_eig$vectors %*% diag(axises) %*% t(a_eig$vectors)


#Sigma
#ggplot() +
#  geom_ellipse(aes(x0 = 0, y0 = 0, a = 4, b = 1, angle = acos(1)),
#               color = 'red4', lty = 2)

#### calculate ellipses #####
#comp_ellipses = function(mean_list, covariance_list, confidence = 0.95){
  
num_elli = dim(mean_list)[1]
x1_list = mean_list[,1]
x2_list = mean_list[,2]
  
a_list = rep(num_elli)
b_list = rep(num_elli)
angle_list = rep(num_elli)

for (j in 1:num_elli) {
  res_eig = eigen(covariance_list[[j]])
  a_list[j] = sqrt(res_eig$values[1])
  b_list[j] = sqrt(res_eig$values[2])
  if(sign(res_eig$vectors[1,1]) == sign(res_eig$vectors[2,2]) ){
    angle_list[j] = acos(res_eig$vectors[1,1])
  }
  else{angle_list[j] = acos(-res_eig$vectors[1,1]) }
}  
  
a_list = a_list * sqrt(qchisq(0.95 ,df=2 ))
b_list = b_list * sqrt(qchisq(0.95 ,df=2 ))
  
########################################### try 1
#library(MVQuickGraphs)
#n = 1
#p = 3
#for (j in 1:num_elli) {
#  confidenceEllipse(X.mean = c(a_list[j], b_list[j] ),
#                    eig = eigen(covariance_list[[j]]),
#                    n = n, p = p,
#                    alpha = 0.05)
#}
########################################### try 1
########################################### try 2
#library(car)
#plot(0,0,xlim=c(-3,0.3),ylim=c(-3,2),  col="red")
#points(x1_list,x2_list)
#for (j in 1:num_elli) {
#  ellipse(center=c(x1_list[j], x2_list[j]),shape=covariance_list[[j]],radius=sqrt(qchisq(0.9,2)), col = "red")
#}
########################################### try 2
##### plot ##### 
library("ggplot2")
library("ggrepel")
library("deldir")
#library("ggvoronoi")
############# plot rows of R ##############################
#plotdata.Rt <- data.frame(R1=Z_final[,1],R2=Z_final[,2])
#p <- ggplot(data = plotdata.Rt, aes(x = R1, y = R2)) + theme_bw() + 
#  geom_point(shape = 16, colour = "gray89", size = 1,stroke = 1) +
#  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#  xlim(-3.5, 3) + ylim(-2.8,2.7) + 
#  xlab("") +
#  ylab("")
#p
############# plot the high-degree authors###############
# plotdata.centers <- data.frame(R1=x1, R2=x2)
# selectAuthors = embeddings$names[embeddings$selected_authors]
# label_point = selected_author_label
# L <- length(selected_author_label)

p1 <- ggplot()+  geom_ellipse(aes(x0 = x1_list, y0 = x2_list, a = a_list, b = b_list, angle =angle_list),
                              color = 'red4', lty = 2) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
# +  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-7, 1.5) + ylim(-6,5) + 
  xlab("") +
  ylab("")
p1

plotdata.select <- data.frame(R1 = x1_list,
                              R2 = x2_list)
set.seed(0)
label_size <- c(rep(4.3,L))
p1 <- p1 + geom_point(data = plotdata.select, aes(x = R1, y = R2),
                    shape = 16,
                    colour = c(rep("blue",nrow(plotdata.select))),
                    size = 1,stroke = 1) + 
  geom_text_repel(data = plotdata.select, 
                  aes(x = R1, y = R2, label = label_point),
                  segment.size = 0.8,
                  colour = c(rep("red",nrow(plotdata.select))),
                  size = label_size, 
                  box.padding = 0.8, max.overlaps = 100)

p1

######################################################
###### Discussion regarding specific authors #########
######################################################

############## confidence region ###############
# num_auth = length(x1_list)
# for (r in 1:num_auth) {
#  p <- p + geom_ellipse(aes(x0 = x1_list[r], y0 = x2_list[r], a = a_list[r], b = b_list[r], angle =angle_list[r]),
#                        color = 'red4', lty = 1)
# }

# p <- p + 
#   geom_ellipse(aes(x0 = x1_list, y0 = x2_list, a = a_list, b = b_list, angle =angle_list),
#                color = 'red4', lty = 1)
# p

# dev.copy2pdf(file=paste("ResearchTriangle.pdf",sep = ''),width = 12,height = 10)


######### Zhao and Liang #############
# which(authors3000 == "Robert Tibshirani")
# which(authors3000 == "Lueping Zhao")
# which(authors3000 == "Kung Yee Liang")
# paper_id1 = which(hypergraph[,5] != 0)
# paper_id2 = which(hypergraph[,56] != 0)
# joint_paper_id = intersect(paper_id1, paper_id2)
# all_paper_id = union(paper_id1, paper_id2)
# embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)
# paper_id1 =setdiff(paper_id1, joint_paper_id)
# paper_id2 =setdiff(paper_id2, joint_paper_id)
# # 
# embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
# embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
# embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
# x1 = embeddings[,1]
# x2 = embeddings[,2]
# Citee = c(rep("Liang, K.Y.", length(paper_id1)),rep("Zhao, L.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
# 
# F_em = final_embeddings$F_em
# Z_em = final_embeddings$Z_em
# 
# Z_em = embeddings$Z_mat
# plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[lky_paper_id,1], F_em[lky_paper_id,2] )
# plot(F_em[zlp_paper_id,1], F_em[zlp_paper_id,2] )
# 
# plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)
# 
# which(selected_author_label == "Lueping Zhao")
# which(selected_author_label == "Kung Yee Liang")
# 
# p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(29,33)], y0 = x2_list[c(29,33)], a = a_list[c(29,33)], 
#                                  b = b_list[c(29,33)], angle =angle_list[c(29,33)]),
#                               color = 'red4', lty = 2) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#   #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlim(-5.6, 1.2) + ylim(-5,4) + 
#   xlab("") +
#   ylab("")
# p
# plotdata.select <- data.frame(x1 = x1_list[c(29,33)],
#                               x2 = x2_list[c(29,33)])
# label_point = selected_author_label[c(29,33)]
# set.seed(0)
# label_size <- c(rep(5,2))
# p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
#                       shape = 16,
#                       colour = c(rep("blue",nrow(plotdata.select))),
#                       size = 1,stroke = 1) + 
#   geom_text_repel(data = plotdata.select, 
#                   aes(x = x1, y = x2, label = label_point),
#                   segment.size = 0.8,
#                   colour = c(rep("red",nrow(plotdata.select))),
#                   size = label_size, 
#                   box.padding = 0.8, max.overlaps = 100)
# 
# p
# 
# 
# 
# cols <- c("#1170AA", "#55AD89", "#EF6F6A")
# 
# p = p + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
#   scale_color_manual(values = cols)
# p
# 

############ Tibshirani and Hastie #############

which(authors3000 == "Robert Tibshirani")
which(authors3000 == "Trevor Hastie")

alpha_est[c(34,48)]
hist(Theta_est[,34])
summary(invlogit(Theta_est[,34]))
summary(invlogit(Theta_est[,2600]))
hist(Theta_est[,48])
summary(invlogit(Theta_est[,48]))
summary(b_hes(Theta_est[,34]))
summary(b_hes(Theta_est[,48]))

paper_id1 = which(hypergraph[,34] != 0)
paper_id2 = which(hypergraph[,48] != 0)
joint_paper_id = intersect(paper_id1, paper_id2)
all_paper_id = union(paper_id1, paper_id2)

embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)

paper_id1 =setdiff(paper_id1, joint_paper_id)
paper_id2 =setdiff(paper_id2, joint_paper_id)

embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
x1 = embeddings[,1]
x2 = embeddings[,2]
upper_list = which(x2>1.3)
all_paper_id[upper_list]
row.names(hypergraph)[all_paper_id[upper_list]]

which(row.names(hypergraph) == "10530") # need to minus 31
hypergraph[2483,34]

citer = PapPapMat$ToPap[which(PapPapMat$FromPap == 17620)]
paper$title[citer]
author_s = AuPapMat$idxAu[which((AuPapMat$idxPap) %in% citer)]

authors_name = read.csv("author_name.txt")
authors_name[author_s,1] 


paper[10530,]

Citee = c(rep("Tibshirani, R.", length(paper_id1)),rep("Hastie, T.", length(paper_id2) ), rep("Both",length(joint_paper_id)))

F_em = final_embeddings$F_em
# Z_em = final_embeddings$Z_em
# Z_em = embeddings$Z_mat
# plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
plot(F_em[paper_id1,1], F_em[paper_id1,2] )
plot(F_em[paper_id2,1], F_em[paper_id2,2] )

plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)

which(selected_author_label == "Robert Tibshirani")
which(selected_author_label == "Trevor Hastie")

p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(15,16)], y0 = x2_list[c(15,16)], a = a_list[c(15,16)], 
                                     b = b_list[c(15,16)], angle =angle_list[c(15,16)]),
                                 color = 'red4', lty = 2) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
  #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-2, -0.5) + ylim(-3,-0.5) + 
  xlab("") +
  ylab("")
p
plotdata.select <- data.frame(x1 = x1_list[c(15,16)],
                              x2 = x2_list[c(15,16)])
label_point = selected_author_label[c(15,16)]
set.seed(0)
label_size <- c(rep(5,2))
p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
                    shape = 16,
                    colour = c(rep("blue",nrow(plotdata.select))),
                    size = 1,stroke = 1) + 
  geom_text_repel(data = plotdata.select, 
                  aes(x = x1, y = x2, label = label_point),
                  segment.size = 0.8,
                  colour = c(rep("red",nrow(plotdata.select))),
                  size = label_size, 
                  box.padding = 0.8, max.overlaps = 100)

p

citer = PapPapMat$ToPap[which(PapPapMat$FromPap == 17620)]
paper$title[citer]

cols <- c("#1170AA", "#55AD89", "#EF6F6A")

ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
  scale_color_manual(values = cols) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("")
# ############ Danyu Y. Lin and Lee-jen Wei #############
# 
# which(authors3000 == "Danyu Y. Lin")
# which(authors3000 == "Lee-jen Wei")
# 
# alpha_est[c(39,20)]
# hist(Theta_est[,39])
# summary(invlogit(Theta_est[,39]))
# hist(Theta_est[,20])
# summary(invlogit(Theta_est[,20]))
# 
# paper_id1 = which(hypergraph[,39] != 0)
# paper_id2 = which(hypergraph[,20] != 0)
# joint_paper_id = intersect(paper_id1, paper_id2)
# all_paper_id = union(paper_id1, paper_id2)
# 
# embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)
# 
# 
# 
# paper_id1 =setdiff(paper_id1, joint_paper_id)
# paper_id2 =setdiff(paper_id2, joint_paper_id)
# 
# embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
# embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
# embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
# x1 = embeddings[,1]
# x2 = embeddings[,2]
# Citee = c(rep("Lin, D.Y.", length(paper_id1)),rep("Wei, L.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
# 
# F_em = final_embeddings$F_em
# # Z_em = final_embeddings$Z_em
# 
# # Z_em = embeddings$Z_mat
# 
# 
# # plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[paper_id1,1], F_em[paper_id1,2] )
# plot(F_em[paper_id2,1], F_em[paper_id2,2] )
# 
# plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)
# 
# which(selected_author_label == "Danyu Y. Lin")
# which(selected_author_label == "Lee-Jen Wei")
# 
# p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(23,22)], y0 = x2_list[c(23,22)], a = a_list[c(23,22)], 
#                                      b = b_list[c(23,22)], angle =angle_list[c(23,22)]),
#                                  color = 'red4', lty = 2) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#   #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlim(-5.6, 1.2) + ylim(-5,4) + 
#   xlab("") +
#   ylab("")
# p
# plotdata.select <- data.frame(x1 = x1_list[c(23,22)],
#                               x2 = x2_list[c(23,22)])
# label_point = selected_author_label[c(23,22)]
# set.seed(0)
# label_size <- c(rep(5,2))
# p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
#                     shape = 16,
#                     colour = c(rep("blue",nrow(plotdata.select))),
#                     size = 1,stroke = 1) + 
#   geom_text_repel(data = plotdata.select, 
#                   aes(x = x1, y = x2, label = label_point),
#                   segment.size = 0.8,
#                   colour = c(rep("red",nrow(plotdata.select))),
#                   size = label_size, 
#                   box.padding = 0.8, max.overlaps = 100)
# 
# p
# 
# 
# 
# cols <- c("#1170AA", "#55AD89", "#EF6F6A")
# 
# ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
#   scale_color_manual(values = cols) + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlab("") +
#   ylab("")
# p
# 
# 
# 
# ############ J.L. Wang and Hastie #############
# 
# which(authors3000 == "Jane-ling Wang")
# which(authors3000 == "Trevor Hastie")
# 
# alpha_est[c(189,48)]
# hist(Theta_est[,189])
# summary(invlogit(Theta_est[,189]))
# hist(Theta_est[,48])
# summary(invlogit(Theta_est[,48]))
# 
# paper_id1 = which(hypergraph[,189] != 0)
# paper_id2 = which(hypergraph[,48] != 0)
# joint_paper_id = intersect(paper_id1, paper_id2)
# all_paper_id = union(paper_id1, paper_id2)
# 
# embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)
# 
# 
# 
# paper_id1 =setdiff(paper_id1, joint_paper_id)
# paper_id2 =setdiff(paper_id2, joint_paper_id)
# 
# embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
# embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
# embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
# x1 = embeddings[,1]
# x2 = embeddings[,2]
# Citee = c(rep("Wang, J.L.", length(paper_id1)),rep("Hastie, T.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
# 
# F_em = final_embeddings$F_em
# # Z_em = final_embeddings$Z_em
# 
# # Z_em = embeddings$Z_mat
# 
# 
# # plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[paper_id1,1], F_em[paper_id1,2] )
# plot(F_em[paper_id2,1], F_em[paper_id2,2] )
# 
# plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)
# 
# which(selected_author_label == "Jane-Ling Wang")
# which(selected_author_label == "Trevor Hastie")
# 
# p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(18,16)], y0 = x2_list[c(18,16)], a = a_list[c(18,16)], 
#                                      b = b_list[c(18,16)], angle =angle_list[c(18,16)]),
#                                  color = 'red4', lty = 2) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#   #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlim(-5.6, 1.2) + ylim(-5,4) + 
#   xlab("") +
#   ylab("")
# p
# plotdata.select <- data.frame(x1 = x1_list[c(18,16)],
#                               x2 = x2_list[c(18,16)])
# label_point = selected_author_label[c(18,16)]
# set.seed(0)
# label_size <- c(rep(5,2))
# p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
#                     shape = 16,
#                     colour = c(rep("blue",nrow(plotdata.select))),
#                     size = 1,stroke = 1) + 
#   geom_text_repel(data = plotdata.select, 
#                   aes(x = x1, y = x2, label = label_point),
#                   segment.size = 0.8,
#                   colour = c(rep("red",nrow(plotdata.select))),
#                   size = label_size, 
#                   box.padding = 0.8, max.overlaps = 100)
# 
# p
# 
# 
# 
# cols <- c("#1170AA", "#55AD89", "#EF6F6A")
# 
# ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
#   scale_color_manual(values = cols) + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlab("") +
#   ylab("")
# p
# 
# 
# 
# 
# 
# ############ Iain Johnstone and Hastie #############
# 
# which(authors3000 == "Iain Johnstone")
# which(authors3000 == "David Donoho")
# 
# alpha_est[c(50,36)]
# hist(Theta_est[,36])
# summary(invlogit(Theta_est[,36]))
# hist(Theta_est[,50])
# summary(invlogit(Theta_est[,50]))
# 
# paper_id1 = which(hypergraph[,50] != 0)
# paper_id2 = which(hypergraph[,36] != 0)
# joint_paper_id = intersect(paper_id1, paper_id2)
# all_paper_id = union(paper_id1, paper_id2)
# 
# embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)
# 
# 
# 
# paper_id1 =setdiff(paper_id1, joint_paper_id)
# paper_id2 =setdiff(paper_id2, joint_paper_id)
# 
# embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
# embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
# embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
# x1 = embeddings[,1]
# x2 = embeddings[,2]
# Citee = c(rep("Joinstone, I.", length(paper_id1)),rep("Donoho, D.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
# 
# F_em = final_embeddings$F_em
# # Z_em = final_embeddings$Z_em
# 
# # Z_em = embeddings$Z_mat
# 
# 
# # plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[paper_id1,1], F_em[paper_id1,2] )
# plot(F_em[paper_id2,1], F_em[paper_id2,2] )
# 
# plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)
# 
# which(selected_author_label == "Iain Johnstone")
# which(selected_author_label == "David Donoho")
# 
# p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(2,1)], y0 = x2_list[c(2,1)], a = a_list[c(2,1)], 
#                                      b = b_list[c(2,1)], angle =angle_list[c(2,1)]),
#                                  color = 'red4', lty = 2) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#   #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlim(-0.15, 0.6) + ylim(-2.5,-1) + 
#   xlab("") +
#   ylab("")
# p
# plotdata.select <- data.frame(x1 = x1_list[c(2,1)],
#                               x2 = x2_list[c(2,1)])
# label_point = selected_author_label[c(2,1)]
# set.seed(0)
# label_size <- c(rep(5,2))
# p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
#                     shape = 16,
#                     colour = c(rep("blue",nrow(plotdata.select))),
#                     size = 1,stroke = 1) + 
#   geom_text_repel(data = plotdata.select, 
#                   aes(x = x1, y = x2, label = label_point),
#                   segment.size = 0.8,
#                   colour = c(rep("red",nrow(plotdata.select))),
#                   size = label_size, 
#                   box.padding = 0.8, max.overlaps = 100)
# 
# p
# 
# 
# 
# cols <- c("#1170AA", "#55AD89", "#EF6F6A")
# 
# ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
#   scale_color_manual(values = cols) + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlab("") +
#   ylab("")
# p
# 
# 
# 
# ############ Iain Johnstone and Hastie #############
# 
# which(authors3000 == "Iain Johnstone")
# which(authors3000 == "David Donoho")
# 
# alpha_est[c(50,36)]
# hist(Theta_est[,36])
# summary(invlogit(Theta_est[,36]))
# hist(Theta_est[,50])
# summary(invlogit(Theta_est[,50]))
# 
# paper_id1 = which(hypergraph[,50] != 0)
# paper_id2 = which(hypergraph[,36] != 0)
# joint_paper_id = intersect(paper_id1, paper_id2)
# all_paper_id = union(paper_id1, paper_id2)
# 
# embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)
# 
# 
# 
# paper_id1 =setdiff(paper_id1, joint_paper_id)
# paper_id2 =setdiff(paper_id2, joint_paper_id)
# 
# embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
# embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
# embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
# x1 = embeddings[,1]
# x2 = embeddings[,2]
# Citee = c(rep("Joinstone, I.", length(paper_id1)),rep("Donoho, D.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
# 
# F_em = final_embeddings$F_em
# # Z_em = final_embeddings$Z_em
# 
# # Z_em = embeddings$Z_mat
# 
# 
# # plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[paper_id1,1], F_em[paper_id1,2] )
# plot(F_em[paper_id2,1], F_em[paper_id2,2] )
# 
# plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)
# 
# which(selected_author_label == "Iain Johnstone")
# which(selected_author_label == "David Donoho")
# 
# p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(2,1)], y0 = x2_list[c(2,1)], a = a_list[c(2,1)], 
#                                      b = b_list[c(2,1)], angle =angle_list[c(2,1)]),
#                                  color = 'red4', lty = 2) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
#   #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlim(-0.15, 0.6) + ylim(-2.5,-1) + 
#   xlab("") +
#   ylab("")
# p
# plotdata.select <- data.frame(x1 = x1_list[c(2,1)],
#                               x2 = x2_list[c(2,1)])
# label_point = selected_author_label[c(2,1)]
# set.seed(0)
# label_size <- c(rep(5,2))
# p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
#                     shape = 16,
#                     colour = c(rep("blue",nrow(plotdata.select))),
#                     size = 1,stroke = 1) + 
#   geom_text_repel(data = plotdata.select, 
#                   aes(x = x1, y = x2, label = label_point),
#                   segment.size = 0.8,
#                   colour = c(rep("red",nrow(plotdata.select))),
#                   size = label_size, 
#                   box.padding = 0.8, max.overlaps = 100)
# 
# p
# 
# 
# 
# cols <- c("#1170AA", "#55AD89", "#EF6F6A")
# 
# ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
#   scale_color_manual(values = cols) + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   xlab("") +
#   ylab("")
# p
# 
# 
# 
# 
# 
# 
# 


############ Luke Tierney and Robert Kass #############

which(authors3000 == "Luke Tierney")
which(authors3000 == "Robert Kass")

alpha_est[c(31,108)]
hist(Theta_est[,31])
summary(invlogit(Theta_est[,108]))
hist(Theta_est[,31])
summary(invlogit(Theta_est[,108]))

paper_id1 = which(hypergraph[,31] != 0)
paper_id2 = which(hypergraph[,108] != 0)
joint_paper_id = intersect(paper_id1, paper_id2)
all_paper_id = union(paper_id1, paper_id2)

embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)

paper_id1 =setdiff(paper_id1, joint_paper_id)
paper_id2 =setdiff(paper_id2, joint_paper_id)

embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
x1 = embeddings[,1]
x2 = embeddings[,2]
Citee = c(rep("Tierney, L.", length(paper_id1)),rep("Kass, R.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
F_em = final_embeddings$F_em
plot(F_em[paper_id1,1], F_em[paper_id1,2] )
plot(F_em[paper_id2,1], F_em[paper_id2,2] )

plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)

which(selected_author_label == "Luke Tierney")
which(selected_author_label == "Robert Kass")

p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(12,19)], y0 = x2_list[c(12,19)], a = a_list[c(12,19)], 
                                     b = b_list[c(12,19)], angle =angle_list[c(12,19)]),
                                 color = 'red4', lty = 2) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
  #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-1.9, -1.32) + ylim(2,2.8) + 
  xlab("") +
  ylab("")
p
plotdata.select <- data.frame(x1 = x1_list[c(12,19)],
                              x2 = x2_list[c(12,19)])
label_point = selected_author_label[c(12,19)]
set.seed(0)
label_size <- c(rep(5,2))
p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
                    shape = 16,
                    colour = c(rep("blue",nrow(plotdata.select))),
                    size = 1,stroke = 1) + 
  geom_text_repel(data = plotdata.select, 
                  aes(x = x1, y = x2, label = label_point),
                  segment.size = 0.8,
                  colour = c(rep("red",nrow(plotdata.select))),
                  size = label_size, 
                  box.padding = 0.8, max.overlaps = 100)

p



cols <- c("#1170AA", "#55AD89", "#EF6F6A")

ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
  scale_color_manual(values = cols) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("")


############ Luis Pericchi and Larry Wasserman #############

which(authors3000 == "Luis Pericchi")
which(authors3000 == "Larry Wasserman")

alpha_est[c(434,190)]
hist(Theta_est[,434])
summary(invlogit(Theta_est[,190]))
hist(Theta_est[,434])
summary(invlogit(Theta_est[,190]))

paper_id1 = which(hypergraph[,434] != 0)
paper_id2 = which(hypergraph[,190] != 0)
joint_paper_id = intersect(paper_id1, paper_id2)
all_paper_id = union(paper_id1, paper_id2)

embeddings = matrix(0, nrow = length(all_paper_id), ncol = 2)

paper_id1 =setdiff(paper_id1, joint_paper_id)
paper_id2 =setdiff(paper_id2, joint_paper_id)

embeddings[1:length(paper_id1),1:2] = F_em[paper_id1,]
embeddings[(length(paper_id1) + 1): (length(paper_id1) + length(paper_id2)),1:2] = F_em[paper_id2,]
embeddings[(length(paper_id1) + length(paper_id2)+ 1):length(all_paper_id) ,1:2] = F_em[joint_paper_id,]
x1 = embeddings[,1]
x2 = embeddings[,2]
Citee = c(rep("Pericchi, L.", length(paper_id1)),rep("Wasserman, L.", length(paper_id2) ), rep("Both",length(joint_paper_id)))
F_em = final_embeddings$F_em
plot(F_em[paper_id1,1], F_em[paper_id1,2] )
plot(F_em[paper_id2,1], F_em[paper_id2,2] )

plot.papers =data.frame(x1 = x1, x2 = x2, Citee = Citee)

which(selected_author_label == "Luis Pericchi")
which(selected_author_label == "Larry Wasserman")

p <- ggplot(NULL)+  geom_ellipse(aes(x0 = x1_list[c(10,20)], y0 = x2_list[c(10,20)], a = a_list[c(10,20)], 
                                     b = b_list[c(10,20)], angle =angle_list[c(10,20)]),
                                 color = 'red4', lty = 2) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #  geom_point(data = plotdata.Rt, aes(x = R1, y = R2),
  #             shape = 16, colour = "gray89", size = 1,stroke = 1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-1.2, 0) + ylim(1.5,3) + 
  xlab("") +
  ylab("")
p
plotdata.select <- data.frame(x1 = x1_list[c(10,20)],
                              x2 = x2_list[c(10,20)])
label_point = selected_author_label[c(10,20)]
set.seed(0)
label_size <- c(rep(5,2))
p <- p + geom_point(data = plotdata.select, aes(x = x1, y = x2),
                    shape = 16,
                    colour = c(rep("blue",nrow(plotdata.select))),
                    size = 1,stroke = 1) + 
  geom_text_repel(data = plotdata.select, 
                  aes(x = x1, y = x2, label = label_point),
                  segment.size = 0.8,
                  colour = c(rep("red",nrow(plotdata.select))),
                  size = label_size, 
                  box.padding = 0.8, max.overlaps = 100)

p



cols <- c("#1170AA", "#55AD89", "#EF6F6A")

ggplot() + geom_point(data=plot.papers, aes(x1, x2, group = Citee,shape=Citee, color=Citee) ) +
  scale_color_manual(values = cols) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("")




