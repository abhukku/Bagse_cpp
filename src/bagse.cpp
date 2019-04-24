#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

using namespace std;
using namespace arma;

//' An Rcpp Armadillo function to find log10 Bayes factors for an individual gene corresponding to a grid of a mixture normal
//'
//' @param bhat          A double that represents an estimated effect size of a single gene
//' @param sdhat         A double that represents the standard error of the estimated effect size of a single gene
//' @param phi_vec       A k length column vector type that represents the standard deviations of the mixture of normal distributions distribution of the true effect size of an associated gene
//' @return              A k length column vector of log10 Bayes Factors, representing the likelihood of the gene coming from that that component of the mixture normal
// [[Rcpp::export]] 
arma::colvec compute_log10BF(double bhat, double sdhat, arma::colvec& phi_vec) {
	int p = phi_vec.size();
	colvec val(p);
	for (int i = 0; i < p;i++) {
		val[i] = 0.5*log10(pow(sdhat,2.0)/(pow(sdhat,2.0)+pow(phi_vec[i],2.0)))+0.5*((pow(bhat/sdhat,2.0)*pow(phi_vec[i],2.0))/(pow(sdhat,2)+pow(phi_vec[i],2.0)))/log(10);
	}
	return val;
}

//' An Rcpp Armadillo function that computes a weighted sum of log10 Bayes Factors 
//' 
//' @param log10_BFv      A k length column vector with log10 Bayes Factors corresponding to components of the mixture normal for the effect size of a gene
//' @param wv             A k length column vector of weights corresponding to the components of the mixture normal for the effect size of a gene
//' @return               The weighted sum of the log10 Bayes Factors of the components of the mixture normal for the effect size of a gene
// [[Rcpp::export]] 
double weightedlog10BF(arma::colvec& log10_BFv, arma::colvec& wv) {
	double max_bf = log10_BFv.max();
	double totsum = 0;
	for (int i = 0; i < wv.size(); i++) {
		totsum += pow(10.0,log10_BFv[i]-max_bf)*wv[i];
	}
	return(log10(totsum)+max_bf);
}
//' Rcpp Armadillo function to compute a vector containing posterior weights for the mixture normal distribution of an associated gene's effect size
//'
//' @param log10_BF_vec           A k length vector of Bayes Factors representing likelihoods of the gene coming from each of the p mixture normal components
//' @param pi_vec                 A k length vector of prior probabilities of the gene coming from each of the p mixture normal components
//' @param pi0                    A double type prior probability of the gene being unassociated 
//' @param log10_NC               A double type number, which is log10 of a normalizing constant, for the logsum of the prior probabilities and the Bayes Factors
//' @return                       A k+1 length vector of posterior probabilities of the gene coming from each component of the mixture normal. The first entry is the posterior probability of the gene being unassociated
// [[Rcpp::export]] 
arma::colvec computegridwts(arma::colvec& log10_BF_vec, arma::colvec& pi_vec, double pi0, double log10_NC) {
	int i,j;
	double sum = 0;
	colvec p_vec(pi_vec.size()+1);
	for (i = 0; i < pi_vec.size(); i++) {
		if (pi_vec[i]<1e-10) {
			pi_vec[i] = 1e-10;
		}
		sum += pi_vec[i];
	}

	colvec lbv(log10_BF_vec.size()+1);
	for (i = 0; i < p_vec.size();i++) {
		if (i == 0) {
			p_vec[i] = pi0;
			lbv[i] = 0;
		}
		else {
			p_vec[i] = (1-pi0)*pi_vec[i-1]/sum;
			lbv[i] = log10_BF_vec[i-1];
		}
	}
	/*
	for (i = 0; i < p_vec.size();i++) {
		if (i == 0) {
			lbv[i] = 0;
		}
		else {
			lbv[i] = log10_BF_vec[i-1];
		}
	}
	*/
	colvec rst(p_vec.size());
	for (j=0;j<rst.size();j++) {
		rst[j] = pow(10.0,log10(p_vec)[j]+lbv[j]-log10_NC);
	}
	return rst;

}
//' An Rcpp Armadillo function to compute an iteration of the key EM algorithm 
//'
//' @param p             A k+2 length column vector, where k is the components in the mixture normal, that contains all the parameters to be estimated, including the enrichment parameters and the mixture normal weights
//' @param BF_matrix     A (n x k) matrix of Bayes Factors, which is the k Bayes Factors representing each mixture normal component for all n genes
//' @param annot_vec     A n length integer vector of 0s and 1s where a 1 represents that the gene is in the gene set of interest
//' @return              A k+2 length column vector with the new parameter estimates, with the first two being the enrichment parameters and the rest being the mixture normal weight estimates
// [[Rcpp::export]] 
arma::colvec torus_pool_em(arma::colvec& p, arma::mat& BF_matrix, arma::uvec& annot_vec) {
	int size = BF_matrix.n_cols;
	int sp = BF_matrix.n_rows;
	uvec cat_vec = unique(annot_vec);
	int enrich_param_size = cat_vec.size();
	vec pi_vec = p.subvec(enrich_param_size,(enrich_param_size+size-1));
	uvec ind(sp);
	colvec log10_BF_vec(sp);
	colvec prior_vec(sp);
	colvec log10_NC_vec(sp);
	vec pip_vec(sp);
	colvec wts(pi_vec.size()+1);
	colvec nwts(pi_vec.size());
	arma::mat wts_matrix(sp,pi_vec.size()+1);
	int i,j; 
	for (i = 0; i < sp; i++) {
		for (j = 0; j < enrich_param_size; j++) {
			if (cat_vec[j] == annot_vec[i]) {
				ind[i] = j;
			}
		}
		colvec bf_vec = conv_to< colvec >::from(BF_matrix.row(i));
		colvec pivec = conv_to< colvec > ::from(pi_vec);
		log10_BF_vec[i] = weightedlog10BF(bf_vec,pivec);
		prior_vec[i] = p[ind[i]];
		colvec wlbffirst(2),wlbfsecond(2);
		wlbffirst[0] = 0; wlbffirst[1] = log10_BF_vec[i];
		wlbfsecond[0] = 1-prior_vec[i]; wlbfsecond[1] = prior_vec[i];
		log10_NC_vec[i] = weightedlog10BF(wlbffirst,wlbfsecond);
		pip_vec[i] = 1-pow(10,log10(1-prior_vec[i])-log10_NC_vec[i]);
		colvec bfrow = conv_to< colvec > ::from(BF_matrix.row(i));

		wts_matrix.row(i) = conv_to< rowvec> ::from(computegridwts(bfrow,pi_vec,1-prior_vec[i],log10_NC_vec[i]));
	}
	colvec v(enrich_param_size);
	for (i = 0; i < enrich_param_size; i++) {
		uvec indices = find(annot_vec == cat_vec[i]);
		v[i] = sum(pip_vec.elem(indices))/indices.size();
	}
	colvec pnew(enrich_param_size+nwts.size());
	for (i = 0; i < wts_matrix.n_cols; i++) {
		wts[i] = mean(wts_matrix.col(i));
		if (i>0) {
			nwts[i-1] = wts[i]/(1-wts[0]);
		}
	}
	for (i = 0; i < enrich_param_size; i++) {
		pnew[i] = v[i];
	}
	for (j = enrich_param_size; j < enrich_param_size+nwts.size();j++) {
		pnew[j] = nwts[j-enrich_param_size];
	}
	return pnew;
}


//' An Rcpp Armadillo function to find the log likelihood using the Bayesian model, with a given set of parameters
//'
//' @param p             A k+2 length column vector, where k is the components in the mixture normal, that contains all the parameters, including the enrichment parameters and the mixture normal weights
//' @param BF_matrix     A (n x k) matrix of Bayes Factors, which is the k Bayes Factors representing each mixture normal component for all n genes
//' @param annot_vec     A n length integer vector of 0s and 1s where a 1 represents that the gene is in the gene set of interest
//' @return              A double type representing the log likelihood with this set of parameters
// [[Rcpp::export]] 
double torus_pool_loglik(arma::vec& p, arma::mat& BF_matrix, arma::uvec& annot_vec) {
	int size = BF_matrix.n_cols;
	int sp = BF_matrix.n_rows;
	uvec cat_vec = unique(annot_vec);
	int enrich_param_size = cat_vec.size();
	vec pi_vec = p.subvec(enrich_param_size,(enrich_param_size+size-1));
	uvec ind(sp);
	colvec log10_BF_vec(sp);
	colvec prior_vec(sp);
	colvec log10_NC_vec(sp);
	double loglik;
	double totsum;
	int i,j; 
	for (i = 0; i < sp; i++) {
		for (j = 0; j < enrich_param_size; j++) {
			if (cat_vec[j] == annot_vec[i]) {
				ind[i] = j;
			}
		}
		colvec bf_vec = conv_to< colvec >::from(BF_matrix.row(i));
		colvec pivec = conv_to< colvec > ::from(pi_vec);
		log10_BF_vec[i] = weightedlog10BF(bf_vec,pivec);
		prior_vec[i] = p[ind[i]];
		colvec wlbffirst(2),wlbfsecond(2);
		wlbffirst[0] = 0; wlbffirst[1] = log10_BF_vec[i];
		wlbfsecond[0] = 1-prior_vec[i]; wlbfsecond[1] = prior_vec[i];
		log10_NC_vec[i] = weightedlog10BF(wlbffirst,wlbfsecond);
	}
	totsum = sum(log10_NC_vec);
	loglik = totsum/(log10(exp(1)));
	return loglik;
}

//' An Rcpp Armadillo function to calculate enrichment parameters, local FDR for all genes, and estimates for the weights of the components for the mixture normal for associated genes
//' @param betahat       An n length column vector containing effect size estimates for all n genes
//' @param sebetahat     An n length column vector containing the standard errors of the effect size estimates for all n genes
//' @param annotation    An n length integer vector of 0s and 1s, with a 1 corresponding to a gene in the gene set of interest
//' @param tol           A double type number corresponding to the threshold for when the EM algorithm should converge. Default is 1e-1
//'
//' @return              A list containing the following variables:
//' * enrichment_est     A 2 length column vector containing the enrichment parameter estimates
//' * effect_est         An column vector containing estimates for the weights of the mixture normal of the effect size for associated genes 
//' * EM                 The raw output from the EM algorithm within this function
//' * ldfr               An n length column vector containing the local FDR for all genes, using enrichment information
// [[Rcpp::export]] 
List torus_cpp(arma::colvec& betahat, arma::colvec& sebetahat, arma::uvec& annotation, double tol = 1e-1) {
	double phi_min = sebetahat.min()/10;
	int i,j;
	int iter = 0;
	int n = betahat.size();
	colvec bse(n);
	//mat BF_matrix(n,n);
	for (i = 0; i < n; i++) {
		bse[i] = pow(betahat[i],2)-pow(sebetahat[i],2);
	}
	double phi_max = 2 * pow(bse.max(),0.5);
	if (phi_max < phi_min) {
		phi_max = 8*phi_min;
	}
	double phi_iter = log(phi_max);
	int size_seq_vec = (log(phi_max)-log(phi_min))/(0.5*log(2)) + 1;
	colvec seq(size_seq_vec);
	while (phi_iter >= log(phi_min)) {
		seq[iter] = exp(phi_iter);
		phi_iter = phi_iter-0.5*log(2);
		iter++;
	}
	seq = sort(seq);
	colvec phi_vec(seq.size()+1);
	phi_vec[0] = seq[0]/pow(2,.5);
	mat BF_matrix(n,phi_vec.size());
	for (i =0; i < seq.size(); i++) {
		phi_vec[i+1] = seq[i]; 
	}
	double phi_size = phi_vec.size();
	colvec pi_vec(phi_size);
	for (i = 0; i < phi_size;i++){
		pi_vec[i] = 1.0/phi_size;
	}
	uvec cat_vec = unique(annotation);
	int enrich_param_size = cat_vec.size();
	colvec p0(enrich_param_size+pi_vec.size());
	for (j = 0; j < (enrich_param_size+pi_vec.size()); j++) {
		if (j < enrich_param_size) {
			p0[j] = 1e-3;
		}
		else {
			p0[j] = pi_vec[j-enrich_param_size];
		}
	}
	rowvec bf_row;
	for (i = 0; i < n; i++) {
		colvec computelog10bf = compute_log10BF(betahat[i], sebetahat[i],phi_vec);
		bf_row = conv_to< rowvec >::from(computelog10bf);
		BF_matrix.row(i) = bf_row;
	}
	Environment SQEM("package:SQUAREM");
	Function sqem = SQEM["squarem"];
	List tolList = List::create(Named("tol")=tol);
	Function torus_pool_em("torus_pool_em");
	Function torus_pool_loglik("torus_pool_loglik");
	List l = sqem(Named("p",p0), Named("BF_matrix",BF_matrix), Named("annot_vec",annotation), Named("fixptfn",torus_pool_em), Named("objfn",torus_pool_loglik), Named("control",tolList));
	colvec p = as<colvec>(l["par"]);
	colvec alpha_vec(enrich_param_size);
	for (i = 0; i < enrich_param_size; i++) {
		alpha_vec[i] = log(p[i]/(1-p[i]));
		if (i > 0)  {
			alpha_vec[i] = alpha_vec[i] - alpha_vec[0];
		}
	}
	int size = BF_matrix.n_cols;
	int sp = BF_matrix.n_rows;
	uvec ind(sp);
	colvec log10_BF_vec(sp);
	colvec prior_vec(sp);
	colvec log10_NC_vec(sp);
	colvec lfdr(sp);
	for (j = 0; j < size; j++) {
		pi_vec[j] = p[enrich_param_size+j];
	}
	for (i = 0; i < sp; i++) {
		for (j = 0; j < enrich_param_size; j++) {
			if (cat_vec[j] == annotation[i]) {
				ind[i] = j;
			}
		}
		prior_vec[i] = p[ind[i]];
		colvec bfvec = conv_to< colvec >::from(BF_matrix.row(i));
		colvec pivec = conv_to< colvec >::from(pi_vec);
		log10_BF_vec[i] = weightedlog10BF(bfvec,pivec);
		colvec wlbffirst(2),wlbfsecond(2);
		wlbffirst[0] = 0; wlbffirst[1] = log10_BF_vec[i];
		wlbfsecond[0] = 1-prior_vec[i]; wlbfsecond[1] = prior_vec[i];
		log10_NC_vec[i] = weightedlog10BF(wlbffirst,wlbfsecond);
		lfdr[i] = pow(10,log10(1-prior_vec[i]) - log10_NC_vec[i]);
	}
	colvec effect_est = pi_vec;
	colvec enrichment_est = alpha_vec;
	return Rcpp::List::create(Rcpp::Named("enrichment_est")=enrichment_est,Rcpp::Named("effect_est") = effect_est, Rcpp::Named("EM") = l, Rcpp::Named("lfdr")=lfdr);
}

