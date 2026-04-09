######################################################
########## construct the hypergraph network ##########
######################################################

# The code below includes some exploratory and experimental sections (commented out).
# To generate the figures in the paper, run only the uncommented sections.

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
############### learn the embeddings #################
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

final_embeddings = list()
final_embeddings$F_em = F_est
final_embeddings$Z_em = Z_est
final_embeddings$alpha = alpha_est

save(final_embeddings, file = "hg1991-2000_3000_final_embeddings.RData" )
