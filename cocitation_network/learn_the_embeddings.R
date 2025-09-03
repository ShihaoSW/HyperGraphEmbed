######################################################
########## construct the hypergraph network ##########
######################################################

rm(list = ls())
load("BibtexInfo.RData")
load("AuthorPaperInfo.RData")
bib_data = load("BibtexInfo.RData")
ap_data = load("AuthorPaperInfo.RData")

load("BibtexInfo.RData", ex <- new.env())
ls.str(ex) 
head(PapPapMat)
head(AuPapMat)
head(paper_author)
Authors = sort(unique(AuPapMat$idxAu))
Papers = sort(unique(AuPapMat$idxPap))

n_total_author = length(unique(AuPapMat$idxAu))

Journals = unique(AuPapMat$journal)

head(paper)
paper_ind = which(paper$year %in% seq(1991,2000,by=1) )
length(unique(paper_ind))

paper_list = list()
empty_citee = c()
for (i in 1:length(paper_ind)) {
  cat(i,"\n")
  citee_paper = paper_ind[i]
  cited_paper = PapPapMat$ToPap[which(PapPapMat$FromPap ==citee_paper)]
  author_list = c()
  if(length(cited_paper) > 0){
    for (j in 1:length(cited_paper)) {
      author_list = union(author_list, paper_author[[cited_paper[j]]]$id)
    }
    author_list = sort(author_list)
  }
  else{
    empty_citee = union(empty_citee,citee_paper)
  }
  paper_list[[i]] = author_list
}

##### find cited author degrees #####
author_cited_count = rep(0, n_total_author)
for (i in 1:length(paper_ind)) {
  cat(i,"\n")
  if(length(paper_list[[i]]) > 0){
    author_cited_count[ paper_list[[i]]] = author_cited_count[ paper_list[[i]]] + 1
  }
}
author_cited_count

na = 3000 # number of authors

length(table(author_cited_count))
final_auth_array = order(-author_cited_count)[1:na]

author_names = read.csv("author_name.txt", header = FALSE)
author_names = c(author_names[,1])
author_names[Authors[final_auth_array]]
author_names[final_auth_array]

authors3000 = author_names[final_auth_array]

##### create the hypergraph matrix #####
final_auth_array
hypergraph = matrix(0, ncol = na, nrow = length(paper_ind))
row.names(hypergraph) = paper_ind
for (i in 1:length(paper_ind)) {
  cat(i,"\n")
  hypergraph[i, which(final_auth_array %in%  paper_list[[i]]) ] = 1
}

##### check the orders #####
edge_orders = rep(0, length(paper_ind))
for (i in 1:length(paper_ind)) {
  edge_orders[i] = sum(hypergraph[i,])
}
hist(edge_orders)
table(edge_orders)

author_orders = rep(0, na)
for (j in 1:na) {
   author_orders[j] = sum(hypergraph[,j])  
}
hist(author_orders)
table(author_orders)

##### remove non-informative edges ##### 
remove_edge = c()
for (i in 1:length(paper_ind)) {
  if(sum(hypergraph[i,]) <= 1) {remove_edge = union(remove_edge, i) }
}
remove_edge
hypergraph = hypergraph[setdiff(1:length(paper_ind), remove_edge), ]

##### check the orders again #####
m = dim(hypergraph)[1]
edge_orders = rep(0, m)
for (i in 1:m) {
  edge_orders[i] = sum(hypergraph[i,])
}
hist(edge_orders, breaks = 100, xlim = c(0,110),
     xlab = "Orders of the hyperlinks", main = "Histogram of hyperlink orders",
     )
table(edge_orders)

author_orders = rep(0, na)
for (j in 1:na) {
  author_orders[j] = sum(hypergraph[,j])  
}
hist(author_orders, breaks = 1000, xlab = "Occurances of vertices",
     main = "Histogtam of occurances of vertices")
table(author_orders)

dim(hypergraph)

save(hypergraph, file = "hg1991-2000_3000.RData" )

######################################################
########## construct the hypergraph network ##########
######################################################

library(RSpectra)
source("functions_limit.R")

load("hg1999-2000_3000.RData")
V = hypergraph
dim(V)

K0 = 100 # look at 100 components first  

Vsp = as(V, "dgCMatrix")
V.svd = svds(Vsp, k = K0)

V.svd$d

plot(V.svd$d[1:50], lty = 1, xlab = "Index of singular value", ylab = "Singular value")

K1 = 32
sigmas = V.svd$d[1:K1]
u = V.svd$u[,1:K1]
v = V.svd$v[,1:K1]
P0 = u %*% diag(sigmas) %*% t(v)
low_p = exp(-10)
up_p = 1 - exp(-10)
P = trim_p(P0, low_p, up_p)
Theta = logit(P)
m = dim(Theta)[1]
n = dim(Theta)[2]
alpha = colMeans(Theta)
Theta = Theta - rep(1,m) %*% t(alpha)


K0 = 100
Theta.svd = svds(Theta, k = K0)
plot(Theta.svd$d[1:100], lty = 1, xlab = "Index of singular value", ylab = "Singular value")
K2 = 2
sigmas = Theta.svd$d[1:K2]
u = Theta.svd$u[,1:K2]
v = Theta.svd$v[,1:K2]


F0 = u %*% diag(sqrt(sigmas))
Z0 = v %*% diag(sqrt(sigmas))
alpha0 = alpha

###############################################################################

result = am_pga(V, F0, Z0, alpha0, nT = 9000)
F_est = result$F_hat
Z_est = result$Z_hat
alpha_est = result$alpha_hat

#### keep updating if needed ####
#### F0 = F_est
#### Z0 = Z_est            
#### alpha0 = alpha_est         
#### result = am_pga(V, F0, Z0, alpha0, nT = 10000)
#### F_est = result$F_hat
#### Z_est = result$Z_hat
#### alpha_est = result$alpha_hat

final_embeddings = list()
final_embeddings$F_em = F_est
final_embeddings$Z_em = Z_est
final_embeddings$alpha = alpha_est

save(final_embeddings, file = "hg1991-2000_3000_final_embeddings.RData" )

# load("hg1991-2000_3000_final_embeddings.RData")


# embeddings = list()
# embeddings$alpha = alpha_est
# embeddings$F_mat = F_est
# embeddings$Z_mat = Z_est
# embeddings$selected_authors = selected_authors_index

# name_authors = authors3000
# final_embeddings$names = name_authors
# top_200_by_count = 1:200
# top_40_by_both = name_authors[union(1:42, which(alpha_est %in% sort(alpha_est, decreasing = TRUE)[1:40] ))]
# top_40_by_both
# embeddings$author_select = union(1:40, which(alpha_est %in% sort(alpha_est, decreasing = TRUE)[1:40] ))

# top_200_by_alpha
# top_five_indices = author_idx[which(author_orders[author_idx] %in% sort(author_orders[author_idx], 
#                                                                        decreasing = TRUE)[1:5] )]

# save(embeddings, file = "hg1991-2000_3000_embeddings.RData" )

# dim(F_est)
# dim(Z_est)

# Theta_est = rep(1,m) %*% t(alpha_est) + F_est %*% t(Z_est)
# Theta_est = Theta_est - rep(1,m) %*% t(alpha_est)

# K0 = 6
# Theta.svd = svds(Theta_est, k = K0)
# plot(Theta.svd$d, lty = 1, xlab = "Index of singular value", ylab = "Singular value")
# K2 = 6
# sigmas = Theta.svd$d[1:K2]
# u = Theta.svd$u[,1:K2]
# v = Theta.svd$v[,1:K2]

# collect the results in a list
# HGres.summary = list()

##### select the authors ##### 
# selectAuthors = matrix(0, nrow = L, ncol = 5)
# for (k in 1:L) {
#   author_idx = which(label_kmeans == k)
#   top_five_indices = author_idx[which(alpha_est %in% sort(alpha_est, decreasing = TRUE)[1:5] )]
#   selectAuthors[k,] = author_names[final_auth_array][top_five_indices]
# }
# HGres.summary$selectAuthors = selectAuthors

# selectAuthors_count = matrix(0, nrow = L, ncol = 5)
# for (k in 1:L) {
#   author_idx = which(label_kmeans == k)
#   top_five_indices = author_idx[which(author_orders[author_idx] %in% sort(author_orders[author_idx], decreasing = TRUE)[1:5] )]
#   selectAuthors_count[k,] = author_names[final_auth_array][top_five_indices[1:5]]
# }
# HGres.summary$selectAuthors_count = selectAuthors_count
# HGres.summary$authorNames = author_names[final_auth_array]

# ldy_paper_id = which(hypergraph[,39] != 0)
# wlj_paper_id = which(hypergraph[,20] != 0)
# joint_paper_id = intersect(ldy_paper_id, wlj_paper_id)
# all_paper_id = 

# ldy_paper_id =setdiff(ldy_paper_id, joint_paper_id)
# wlj_paper_id =setdiff(wlj_paper_id, joint_paper_id)

# F_em = final_embeddings$F_em
# Z_em = final_embeddings$Z_em

# Z_em = embeddings$Z_mat

# plot(Z_em[c(20,39),1], Z_em[c(20,39),2] )
# plot(F_em[ldy_paper_id,1], F_em[ldy_paper_id,2] )
# plot(F_em[wlj_paper_id,1], F_em[wlj_paper_id,2] )

