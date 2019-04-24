# BAGSE_cpp: An Rcpp Armadillo Implementation of BAGSE
This repository contains an Rcpp Armadillo Implementation of a Bayesian Analysis of Gene Set Enrichment (BAGSE). BAGSE outputs enrichment parameters, which quantify the level of enrichment of a gene set of interest. 

Use the torus_cpp function within the bagse.cpp file to obtain gene set enrichment parameter estimates. You must input annotations (i.e. a vector of 0s and 1s to mark genes in your gene set of interest) and summary statistics (i.e. z-scores or estimated beta values and their standard errors) from a study such as GWAS or differential expression analyses. 

Note: Must have the packages Rcpp, RcppArmadillo, and SQUAREM for use. 

Refer to the Examples folder with a case using simulated data for more details and clarity on the utilization of this software. 

Original implementations of BAGSE in R and C located at https://github.com/xqwen/bagse. Refer there for more details 

Refer any questions to abhukku@umich.edu
