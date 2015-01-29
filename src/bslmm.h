/*
 Genome-wide Efficient Mixed Model Association (GEMMA)
 Copyright (C) 2011  Xiang Zhou
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __BSLMM_H__                
#define __BSLMM_H__

#include <vector>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif


using namespace std;





class BSLMM {

public:
    
   // model structure
    struct LModel{
        
        HYPBSLMM cHyp;
        vector<size_t> rank;
        
        gsl_matrix * Xgamma;
        gsl_matrix * XtX;
        gsl_vector * Xty;
        gsl_vector * beta;
        gsl_vector * Xbeta;
        
        gsl_matrix_view Xgamma_sub;
        gsl_matrix_view XtX_sub;
        gsl_vector_view Xty_sub;
        gsl_vector_view beta_sub;
        
        void InitialVar(size_t &ni_test, size_t &s_max);
        void AssignVar(uchar **X, gsl_vector *y, map<size_t, size_t> &mapRank2pos, size_t &ns_test);
        
        void getSubVar(size_t size_s);
        double CalcPosterior(const double yty, size_t &ni_test, size_t &ns_test, gsl_rng *gsl_r, int &a_mode, const double &trace_G);
        
        void SetXgamma (uchar **X, map<size_t, size_t> &mapRank2pos, size_t &ns_test);
        void Copy(LModel &model);
        void FreeMem();
        
    };
    
    //multiple function related parameters
    size_t n_type;
    vector<size_t> mFunc; // # of variants of each variant type
    double e, e_shape, e_rate; //hyper parameter in the prior gamma distribution
    double vscale;
    map<string, int> mapFunc2Code;
    int iniType;
    vector <double> theta;
    vector <double> subvar;

    
	// IO related parameters
    size_t UnCompBufferSize;
    vector <size_t> CompBuffSizeVec;
	int a_mode;
	size_t d_pace;
	
	string file_bfile;
	string file_geno;
    string file_vcf;
	string file_out;
	
	// LMM related parameters
	double l_min;
	double l_max;
	size_t n_region;
	double pve_null;
	double pheno_mean;
	
	// BSLMM MCMC related parameters
    
    //JY added win, Wvar, ns_neib;
    size_t win, ns_neib;
    size_t nadd_accept, ndel_accept, nswitch_accept, nother_accept;
    size_t nadd, ndel, nswitch, nother;
    int Switch_Flag;
    
	double h_min, h_max, h_scale;			//priors for h
	double rho_min, rho_max, rho_scale;		//priors for rho
	double logp_min, logp_max, logp_scale;		//priors for log(pi)
	size_t s_min, s_max;			//minimum and maximum number of gammas
	size_t w_step;					//number of warm up/burn in iterations
	size_t s_step;					//number of sampling iterations
	size_t r_pace;					//record pace
	size_t w_pace;					//write pace
	size_t n_accept;				//number of acceptance
	size_t n_mh;					//number of MH steps within each iteration
	double geo_mean;				//mean of the geometric distribution
	long int randseed;
	double trace_G;	
	
	HYPBSLMM cHyp_initial;

	// Summary statistics
	size_t ni_total, ns_total;	//number of total individuals and snps
	size_t ni_test, ns_test;	//number of individuals and snps used for analysis
	size_t n_cvt;				//number of covariates
	double time_UtZ;
	double time_Omega;		//time spent on optimization iterations
	double time_Proposal;        //time spent on constructing the proposal distribution for gamma (i.e. lmm or lm analysis)
	vector<bool> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
	vector<bool> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	
	vector<SNPINFO> snpInfo;		//record SNP information
	
	// Not included in PARAM
	gsl_rng *gsl_r;
	gsl_ran_discrete_t *gsl_t; //JY added dynamic gsl_s
    
	map<size_t, size_t> mapRank2pos;
	map<size_t, size_t> mapOrder2pos; // JY: map order index to snp position
    map<size_t, size_t> mapPos2Order; // JY: map position to snp order
    map<size_t, size_t> mapPos2Rank; // JY: map position to snp rank
    map<size_t, size_t> mapRank2Order; // JY: map rank index to snp order
    map<size_t, size_t> mapOrder2Rank; // JY: map order index to snp rank

	
	// Main Functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	
	void RidgeR(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *eval, const double lambda);

    
	void WriteLog ();
	void WriteLR ();
	void WriteBV (const gsl_vector *bv);
	void WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr);
	void WriteParam (const gsl_vector *alpha);
	void WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const gsl_matrix *Result_theta, const gsl_matrix *Result_sigma, const size_t w_col);
	
	//Subfunctions inside MCMC
	void CalcPgamma (double *p_gamma);
    
	double CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2);
    
	//void InitialMCMC (uchar **X, const gsl_vector *Uty, LModel &model_old, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos);
    
	double CalcPosterior (const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *Utu, gsl_vector *alpha_prime, class HYPBSLMM &cHyp);
	double CalcPosterior (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *UtXb, gsl_vector *Utu, gsl_vector *alpha_prime, gsl_vector *beta, class HYPBSLMM &cHyp);
	void CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp);
	void CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *UtXb, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp);
	double CalcREMLE (const gsl_matrix *Utw, const gsl_vector *Uty, const gsl_vector *K_eval);
	double CalcLR (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, vector<pair<size_t, double> > &loglr_sort);		//calculate the maximum marginal likelihood ratio for each analyzed SNPs with gemma, use it to rank SNPs
	void SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z);
	double ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
	double ProposePi (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
    

	void CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty);
    
	//double CalcPosterior (const double yty, class HYPBSLMM &cHyp);
	//double CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const size_t s_size, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp);
    
	void CalcCC_PVEnZ (gsl_vector *z_hat, class HYPBSLMM &cHyp);
	void CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp);
    
    //JY added
    void WriteIniRank (const vector<string> &iniRank);

    void CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX_gamma, const gsl_vector *Xtz_gamma, gsl_vector *z_res, const size_t &s_size, const double &ztz);
    
    double CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec);
    
    gsl_ran_discrete_t * MakeProposal(const size_t &o, double *p_BF, uchar **X, const gsl_vector *z_res, const map<size_t, int> &mapRank2in);

    //JY revised
    double ProposeGamma (LModel &model_old, LModel &rank_new, const double *p_gamma, const size_t &repeat, uchar **X,  gsl_vector *z, const double &ztz, int &flag_gamma);
    

    void NormRes(gsl_vector * z_res);
    
    //MCMC sub-functions
    //void MCMC (uchar **X_Genotype, gsl_vector *z);

    void MCMC_Free_WorkVar(gsl_matrix *Result_hyp, gsl_matrix *Result_gamma, gsl_vector *z_hat);
    
    void CreateSnpPosVec(vector<SNPPOS> &snp_pos);
    
    void InitialMap(const vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos);
    
    void CreateGammaProposal(const double *p_gamma);
    
    void MHmove(const bool &accept, const int &flag_gamma, LModel &model_old, LModel &model_new, double &logPost_old, double &logPost_new);
    
    void MHsave(const size_t &t, size_t &w, gsl_matrix *Result_hyp, gsl_matrix *Result_gamma, LModel &model_old, vector<pair<double, double> > &beta_g, const double &mean_z, const vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos, double &logPost_old);
    
    double MHPropose(uchar **X_Genotype, const double *p_gamma, gsl_vector *z, const double &ztz, LModel &model_old, LModel &model_new, int &flag_gamma, double &logPost_new, double &logPost_old);
    
    //propose gamma sub-functions
    void AddMarker(double &logp, map<size_t, int> &mapRank2in, class HYPBSLMM &cHyp_new, vector<size_t> &rank_new, const double *p_gamma);
    
    void DelMarker(double &logp, map<size_t, int> &mapRank2in, class HYPBSLMM &cHyp_new, vector<size_t> &rank_new, const double *p_gamma);
    
    void SwitchMarker(double &logp, map<size_t, int> &mapRank2in, LModel &model_old, LModel &model_new, uchar **X, gsl_vector *z, const double &ztz);
    
    void AssignRank(map<size_t, int> &mapRank2in, vector<size_t> &rank_new, const vector<size_t> &rank_old, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new);
    
    void SetXgamma (gsl_matrix *Xgamma, uchar **X, vector<size_t> &rank);
    
    void SetXgamma (LModel &model_old, LModel &model_new, uchar **X, const gsl_vector *y);
    
    void WriteMatrix(const gsl_matrix * X, const string &filename);
    void WriteVector(const gsl_vector * X, const string &filename);
    
	//utility functions
//	double vec_sum (gsl_vector *v);
//	void vec_center (gsl_vector *v);
//	double calc_var (gsl_vector *v);
//	void calc_sigma (MCMC &cMcmc);
//	bool comp_lr (pair<size_t, double> a, pair<size_t, double> b);
    
    void InitialMCMC ( uchar **UtX, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos);
    void SetXgamma ( uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, const vector<size_t> &rank_new, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new);
    double CalcPosterior (const double yty, class HYPBSLMM &cHyp, const gsl_vector *pi_vec, const vector<size_t> &rank);
    double CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *pi_vec, gsl_vector *sigma_vec, const vector<size_t> &rank);

    double ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat,  uchar **X, const gsl_vector *z, const gsl_matrix *Xgamma_old, const gsl_matrix *XtX_old, const gsl_vector *Xtz_old, const double &ztz, int &flag_gamma);
    void WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr);
    void WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w);
    void MCMC (uchar **X, const gsl_vector *y, bool original_method);
    
    // added function for newmodel
    void CalcPivec(const vector<double> &theta, gsl_vector *pi_vec, const vector<SNPPOS> &snp_pos);
    void CalcSvec(const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos);
    double ProposeTheta (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
    double ProposeSubvar (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);

    double ProposeTheta (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, size_t j, gsl_vector *pi_vec_old, gsl_vector *pi_vec_new, const vector<SNPPOS> &snp_pos);
    double ProposeSubvar (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, size_t j, gsl_vector *sigma_vec_new, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos);
    
    double CalcLikelihood (const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, class HYPBSLMM &cHyp, gsl_vector *sigma_vec, const vector<size_t> &rank);
    double CalcPsubvar (const class HYPBSLMM &cHyp);
    double CalcPtheta (const class HYPBSLMM &cHyp);
    double CalcPsubvar (const class HYPBSLMM &cHyp, size_t j);
    double CalcPtheta (const class HYPBSLMM &cHyp, size_t j);
    double CalcLikegamma(const gsl_vector *pi_vec, const vector<size_t> &rank);
    
    double CalcSigma(const class HYPBSLMM &cHyp, const size_t &order_i, const vector<SNPPOS> &snp_pos);
    
};

void PrintVector(const gsl_vector * x);
void PrintVector(const vector <double> &x);
void PrintVector(const vector <size_t> &x);
void PrintVector(const double *x);
void PrintVector(const uchar *x, const size_t length);
void PrintVector(const gsl_vector * x, const size_t s);
void PrintMatrix(const gsl_matrix * X, const size_t nrow, const size_t ncol);

void expVector(vector<double> &expvec, vector<double> &logvec);
void CalcXVbeta(gsl_matrix *X, const gsl_vector * sigma_vec);
#endif


