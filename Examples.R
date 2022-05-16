################################################
## install packages and load into R first
################################################
## install.packages("R2OpenBUGS")
## install.packages("MASS")
## install.packages("coda")
## install.packages("gdata")
## install.packages("igraph")
################################################
library(R2OpenBUGS)
library(MASS)
library(coda)
library(gdata)
library(igraph)
################################################

################################################
## Preload PRIDE_function.R
################################################
source("PRIDE_Function.R")
################################################
## Preload Performance_Criteria.R
################################################
source("Performance_Criteria.R")
################################################

################################################
## Example
## AR(1) vs. AR(2)
## set p = 5
## set N1 = N2 = 50
################################################
################################################
## first group
## AR(2) network
## number of edges: 2*(p-2) + 1
## precision matrix setting
## first order: 0.3
## second order: 0.22
## others: 0
################################################

## number of genes
N.P <- 5

## number of possible interactions
Total.Edges <- choose(N.P, 2)

## set precision matrix
AR2.M <- matrix(0, N.P, N.P)

for(i in 1:(N.P-2)){

AR2.M[i, i+1] <- 0.3
AR2.M[i, i+2] <- 0.22

}

AR2.M[(N.P-1), N.P] <- 0.3

AR2.M <- AR2.M + t(AR2.M)

diag(AR2.M) <- 1

isSymmetric(AR2.M)

COR.M1 <- cov2cor(solve(AR2.M))

Adj.1 <- apply(AR2.M, 1, function(input) ifelse(input == 0, 0, 1))

diag(Adj.1) <- 0

## graph for AR(2) network
network1 <- graph_from_adjacency_matrix(Adj.1, mode = 'undirected', diag = F)

plot(network1, layout = layout.circle)

################################################
## Second group
## AR(1) network
## number of edges = (p-1)
## precision matrix setting
## first order: 0.3
## others: 0
################################################

## set precision matrix
AR1.M <- matrix(0, N.P, N.P)

for(i in 1:(N.P-2)){

AR1.M[i, i+1] <- 0.3

}

AR1.M[(N.P-1), N.P] <- 0.3

AR1.M <- AR1.M + t(AR1.M)

diag(AR1.M) <- 1

isSymmetric(AR1.M)

COR.M2 <- cov2cor(solve(AR1.M))

Adj.2 <- apply(AR1.M, 1, function(input) ifelse(input == 0, 0, 1))

diag(Adj.2) <- 0

## graph for AR(1) network
network2 <- graph_from_adjacency_matrix(Adj.2, mode = 'undirected', diag = F)

plot(network2, layout = layout.circle)

################################################
## Differential network
## number of differential edges = (p-2)
################################################

Diff <- Adj.1 - Adj.2

N.Diff.edge <- sum(upperTriangle(Diff))

Diff.Net <- graph_from_adjacency_matrix(Diff, mode = 'undirected', diag = F)

## graph for differential network
plot(Diff.Net, layout = layout.circle)

## True differential edges
True.Diff <- upperTriangle(Diff, byrow = T)

################################################
## generate data from MVN distribution
################################################

## number of samples in each group
N.S <- 50

## response variable
## Y = 1: case group (AR(2))
## Y = 0: control group (AR(1))
Y <- rep(c(1,0), each = N.S)

Data.1 <- mvrnorm(N.S, rep(0, N.P), COR.M1)

Data.2 <- mvrnorm(N.S, rep(0, N.P), COR.M2)

## covariate matrix
Data.X <- rbind(Data.1, Data.2)


################################################
## no screening version
## estimate all interactions
################################################

Edge.S <- rep(1, choose(N.P, 2))

################################################
## run PRIDE model
################################################

star <- Sys.time()

Output <- PRIDE(N1 = N.S, N2 = N.S, Data.X, Y, Edge.S)

as.numeric(Sys.time()-star)

################################################
## outputs
################################################

## posterior samples
Post.Samp <- Output[[1]]

## generate 5000 posterior samples
## 10 BETA + 10 Gamma + 1 deviance
dim(Post.Samp)

colnames(Post.Samp)

## posterior mean of each variable
summary(Output)

## estimated prob(diff edge)
Prob.diff.edge <- summary(Output)[[1]][(Total.Edges+1):(2*Total.Edges), 1]

## posterior mean of BETA
Est.Inter <- summary(Output)[[1]][1:Total.Edges, 1]

## select differential edges 
## prob(diff edge) >= 0.5 or not
Est.diff.edge <- ifelse(Prob.diff.edge >= 0.5, 1, 0)

Est.Diff.Net <- matrix(0, N.P, N.P)

upperTriangle(Est.Diff.Net, byrow = T) <- Est.diff.edge

Est.Diff.Net <- Est.Diff.Net + t(Est.Diff.Net)

Est.G <- graph_from_adjacency_matrix(Est.Diff.Net, mode = "undirected", diag = F)

## visualize the estimated differential network
plot(Est.G)

################################################
## performance
################################################

ACC.Criterion(rbind(Est.diff.edge, True.Diff))


##############################################################
################################################
## informative screening version
## only estimates those true differential edges
################################################
##############################################################

Edge.S <- True.Diff

N.Est <- sum(Edge.S)

################################################
## run PRIDE model
################################################

star <- Sys.time()

Output <- PRIDE(N1 = N.S, N2 = N.S, Data.X, Y, Edge.S)

as.numeric(Sys.time()-star)

################################################
## outputs
################################################

## posterior samples
Post.Samp <- Output[[1]]

dim(Post.Samp)

colnames(Post.Samp)

## posterior mean of each variable
summary(Output)

## estimated prob(diff edge)
Prob.diff.edge <- summary(Output)[[1]][(N.Est+1):(2*N.Est), 1]

## select differential edges 
## prob(diff edge) >= 0.5 or not
Est.diff.edge <- rep(0, choose(N.P, 2))

Est.diff.edge[which(Edge.S == 1)] <- ifelse(Prob.diff.edge >= 0.5, 1, 0)

Est.Diff.Net <- matrix(0, N.P, N.P)

upperTriangle(Est.Diff.Net, byrow = T) <- Est.diff.edge

Est.Diff.Net <- Est.Diff.Net + t(Est.Diff.Net)

Est.G <- graph_from_adjacency_matrix(Est.Diff.Net, mode = "undirected", diag = F)

## visualize the estimated differential network
plot(Est.G)

################################################
## performance
################################################

ACC.Criterion(rbind(Est.diff.edge, True.Diff))





