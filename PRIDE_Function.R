##########################################################################################################
## PRobabilistic Interactions for Differential Edges (PRIDE)
## function: PRIDE()
## PRIDE is a function to generate probabilistic estimation for differential edges
##########################################################################################################
## Inputs for PRIDE function
## 1. N1; 2. N2; 3.Data.X; 4. Y; 5. Edge.S
##########################################################################################################
## 1. N1 represents the # of samples from the first group 
## 2. N2 represents the # of samples from the second group 
##########################################################################################################
## 3. Data.X is a (N1+N2) by p matrix
## p represents the # of variables (gene here)
## Data.X is the data matrix storing covariates (no response variable)
## Note: need to combine two group covariate data (row bind) before inputting to PRIDE function
##########################################################################################################
## 4. Y: a (N1+N2) vector 
## the element in Y is either 1 (for the case group) or 0 (for the control group)
##########################################################################################################
## 5. Edge.S: a p*(p-1)/2 vector which represents the information from the screening step
## each element in Edge.S is a binary value (0/1)
## Edge.S[i] = 1: i_th interaction needed to be estimated in PRIDE
## Edge.S[i] = 0: i_th interaction do not need to estimate in PRIDE
## If all elements in Edge.S = 1, it means that we need to estimate all p*(p-1)/2 interactions in PRIDE
##########################################################################################################
##########################################################################################################
## The output of PRIDE function
##########################################################################################################
## Because we use the "coda" package to store the posterior samples, 
## the output of the PRIDE function is a "MCMC list" object.
## It contains "T" posterior samples for each parameter in the generating list object.
## "T" can be set by users.
## Users can also define which parameters they are interested in making an inference.
## The default setting is generating 5000 posterior samples for both "BETA" and "Gamma" in PRIDE model
##########################################################################################################

PRIDE <- function(N1, N2, Data.X, Y, Edge.S){

#####################################################
## Data.X
## gene expression data
## (N1+N2)*p matrix
## row represents samples
## col represents genes
## combine case and control
## number of row = (N1+N2)
#####################################################
#####################################################
## standardize each group separately
## standerdize by col (gene)
#####################################################

Data1 <- apply(Data.X[1:N1,], 2, function(input) (input-mean(input))/sd(input))

Data2 <- apply(Data.X[(N1+1):(N1+N2),], 2, function(input) (input-mean(input))/sd(input))

X <- as.matrix(rbind(Data1, Data2))

## number of total sample size used in the analysis
Num <- N1+N2

## number of genes
P <- dim(X)[2]

## number of possible interaction terms
N.P <- choose(P, 2)

## number of interactions to estimate in PRIDE
S.P <- sum(Edge.S)

#####################################################
## construct data matrix for storing interaction terms
#####################################################

XTX <- matrix(0, Num, N.P)

for (i in 1:Num){

XTX[i,] <- upperTriangle(X[i,]%*%t(X[i,]), byrow = T)

}

#####################################################
## select which interactions needed to be estimated
## XTX is the input covariate matrix
#####################################################

XTX <- XTX[, which(Edge.S == 1)]

#####################################################
## Start OpenBUGS code
#####################################################

Bayes.Model <- function() {

#####################################################
## Likelihood
## logistic regression: intercept + interactions
#####################################################

for (i in 1:Num) {

logit(p[i]) <- Alpha + inprod(XTX[i,], BETA[1:S.P])

Y[i] ~ dbern(p[i])

}

#####################################################
## prior distribution for intercept
#####################################################

Alpha ~ dnorm(0, 0.001)

#####################################################
## Spike-and-Slab Lasso prior for BETA
## slab distribution: DE(0, 2)
## spike distribution: DE(0, 20)
## Bernoulli prior distribution on Gamma
## set hyperprior p in Bernoulli distribution = 0.7
#####################################################

for (j in 1:S.P){

Gamma[j] ~ dbern(0.7)

tau[j] <- (18*(equals(Gamma[j], 0))) + 2

BETA[j] ~ ddexp(0, tau[j])

}

}

#####################################################
## set the working directory
## set the model.file used in OpenBUGS
#####################################################

Direction <- getwd()

OpenBUGS.file <- gsub(" ", "", paste(Direction,"/Bayes.Model.odc"))

write.model(Bayes.Model, con = OpenBUGS.file)

my.data <- list("Y", "XTX", "Num", "S.P")

inits <- function() {
list(BETA = rep(0, S.P), Gamma = rep(0, S.P))
}

params <- c("BETA", "Gamma")

#########################################################
## start sampler
## for each parameter, generated 5000 posterior samples
## thin = 10
## burn in = 5000
## iteration = 10000
## users can change the settings here!
#########################################################

out <- bugs(data = my.data, inits = inits, parameters.to.save = params, 
model.file = OpenBUGS.file, codaPkg = TRUE,
n.iter = 10000, n.chains = 1, n.burnin = 5000, n.thin = 10, debug = F)

out.coda <- read.bugs(out)

return(out.coda)

}








