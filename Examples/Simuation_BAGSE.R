### Example for estimating the enrichment parameters using bagse_cpp, with simulated data. 
### The file with the simulated data is "Simulated_Data.txt", in this folder


library(Rcpp)
library(RcppArmadillo)
library(SQUAREM)
sourceCpp("bagse.cpp") #The .cpp file in the src folder here

#10,000 genes
n = 10000

#20% of genes are in gene set of interest
prop_annot = .2

#From real data, this is a realistic estimate of alpha_0
a0 = -2.5

#Simulate this data set for an enrichment parameter alpha_1 of 1
a1 = 1
set.seed(815)

#Generate vector of annotations
annot_vec = c(rep(0,(1-prop_annot)*n),rep(1,prop_annot*n))

#Probabilities for each gene that they are associated, based on enrichment parameters
pi = exp(a0+a1*annot_vec)/(1+exp(a0+a1*annot_vec))

#Generate true association status of all genes, based on pi
gamma = rbinom(n,1,pi)

#Generate z-scores for all genes, depending on if they are truly associated
zscores = rnorm(n)
amp = 4.5
#Based on real data, a long-tailed distribution like t can be representatitve of the distribution of z-scores of associated genes
zscores[gamma == 1] = amp*rt(sum(gamma),df = 10)

#Set the parameters to be used in torus_cpp
betahat = zscores
sebetahat = rep(1,n)
annotation = annot_vec

t = torus_cpp(betahat,sebetahat,annotation)
t$enrichment_est
#Should be -3.06 for a0 and 0.93 for a1. 

#If you want to compare the results to the R package, the R code is in the R folder in the repo
source("bagse.R")
r = torus(betahat,sebetahat,annotation)
r$enrichment_est
#Should be -3.07 for a0 and 0.93 for a1 again, but markedly slower



