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

#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include <limits>


#include "lapack.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "bslmm_float.h"
#include "lmm_float.h"  //for class FUNC_PARAM and MatrixCalcLR
#include "lm_float.h"
#include "mathfunc_float.h"  //for function CenterVector
#else
#include "param.h"
#include "bslmm.h"
#include "lmm.h"
#include "lm.h"
#include "mathfunc.h"
#endif

using namespace std;




void BSLMM::CopyFromParam (PARAM &cPar) 
{
    n_type = n_type;
    mFunc = cPar.mFunc;
    
    UnCompBufferSize = cPar.UnCompBufferSize;
    CompBuffSizeVec = cPar.CompBuffSizeVec;
    win = cPar.win;
    nadd_accept = cPar.nadd_accept;
    ndel_accept= cPar.ndel_accept;
    nswitch_accept= cPar.nswitch_accept;
    nother_accept= cPar.nother_accept;
    nadd= cPar.nadd;
    ndel= cPar.ndel;
    nswitch= cPar.nswitch;
    nother= cPar.nother;
    
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
    file_vcf = cPar.file_vcf;
	file_out=cPar.file_out;
	
	l_min=cPar.h_min;	
	l_max=cPar.h_max;  
	n_region=cPar.n_region;	
	pve_null=cPar.pve_null;
	pheno_mean=cPar.pheno_mean;
	
	time_UtZ=0.0;
	time_Omega=0.0;
	n_accept=0;
    Switch_Flag = 0;
	
	h_min=cPar.h_min;	
	h_max=cPar.h_max;  
	h_scale=cPar.h_scale;
	rho_min=cPar.rho_min;	
	rho_max=cPar.rho_max;  
	rho_scale=cPar.rho_scale;
	logp_min=cPar.logp_min;	
	logp_max=cPar.logp_max;  
	logp_scale=cPar.logp_scale;
	
	s_min=cPar.s_min;
	s_max=cPar.s_max;
	w_step=cPar.w_step;
	s_step=cPar.s_step;
	r_pace=cPar.r_pace;
	w_pace=cPar.w_pace;
	n_mh=cPar.n_mh;
	geo_mean=cPar.geo_mean;
	randseed=cPar.randseed;
	trace_G=cPar.trace_G;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;
	
	return;
}


void BSLMM::CopyToParam (PARAM &cPar) 
{
	cPar.time_UtZ=time_UtZ;
	cPar.time_Omega=time_Omega;
	cPar.time_Proposal=time_Proposal;
	cPar.cHyp_initial=cHyp_initial;
	cPar.n_accept=n_accept;
	cPar.pheno_mean=pheno_mean;
	cPar.randseed=randseed;
    
    cPar.nadd_accept=nadd_accept;
    cPar.ndel_accept=ndel_accept;
    cPar.nswitch_accept=nswitch_accept;
    cPar.nother_accept=nother_accept;
    cPar.nadd=nadd;
    cPar.ndel=ndel;
    cPar.nswitch=nswitch;
    cPar.nother=nother;
	
	return;
}

//JY add to write initial significant SNP id out
void BSLMM::WriteIniRank (const vector<string> &iniRank)
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".inirank.txt";
    
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	for (size_t i=0; i<iniRank.size(); ++i) {

			outfile<<iniRank[i]<<endl;

	}
	
	outfile.clear();
	outfile.close();
	return;
}

void BSLMM::WriteMatrix(const gsl_matrix * X, const string &filename){
    string file_str = "./output/"+file_out;
    file_str += filename;
    ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
    for(size_t i=0; i<X->size1; ++i){
        for(size_t j = 0; j < X->size2; ++j){
            outfile << scientific << setprecision(6) << gsl_matrix_get(X, i, j) << "\t" ;
        }
        outfile << endl;
    }
    
    outfile.clear();
	outfile.close();
    return;
} //write gsl_matrix X with filename = ***.txt

void BSLMM::WriteVector(const gsl_vector * X, const string &filename){
    string file_str = "./output/"+file_out;
    file_str += filename;
    ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
    for(size_t i=0; i<X->size; ++i){
        outfile << scientific << setprecision(6) << gsl_vector_get(X, i)<< endl;
    }
    
    outfile.clear();
	outfile.close();
    return;
}//write gsl_vector X with filename = ***.txt


void PrintVector(gsl_vector * x){
    for(size_t i=0; i < x->size; ++i){
        cout << gsl_vector_get(x, i) << ", ";
    }
    cout << endl; 
}

void PrintVector(gsl_vector * x, size_t s){
    for(size_t i=0; i < s; ++i){
        cout << gsl_vector_get(x, i) << ", ";
    }
    cout << endl;
}

void PrintMatrix(gsl_matrix * X, size_t nrow, size_t ncol){
    for (size_t i=0; i<nrow; i++) {
        gsl_vector_view row = gsl_matrix_subrow(X, i, 0, ncol);
        PrintVector(&row.vector);
    }
}

void PrintVector(vector <double> &x){
    for(size_t i=0; i<x.size(); ++i){
        cout << x[i] << ", ";
    }
    cout << endl; 
}

void PrintVector(vector <size_t> &x){
    for(size_t i=0; i<x.size(); ++i){
        cout << x[i] << ", ";
    }
    cout << endl;
}

void PrintVector(double *x){
    for(size_t i=0; i<41; ++i){
        cout << x[i] << ", ";
    }
    cout << endl; 
}

void PrintVector(uchar *x, size_t length){
    for(size_t i=0; i<length; ++i){
        cout << (int)x[i] << ", ";
    }
    cout << endl;
}

// JY add function to write function

void BSLMM::WriteBV (const gsl_vector *bv) 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".bv.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	size_t t=0;
	for (size_t i=0; i<ni_total; ++i) {
		if (indicator_idv[i]==0) {
			outfile<<"NA"<<endl;
		}		
		else {
			outfile<<scientific<<setprecision(6)<<gsl_vector_get(bv, t)<<endl;
			t++;
		}
	}		
	
	outfile.clear();	
	outfile.close();
	return;
}


void BSLMM::WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w, const vector<snpPos> &snp_pos, const vector<pair<size_t, double> > &pos_loglr)
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"geno_pos"<<"\t"<<"markerID"<<"\t"<<"chr"<<"\t"
			<<"bp"<<"\t"<<"alpha"<<"\t" << "lrt" << "\t"
			<<"beta"<<"\t"<<"gamma" << endl; //JY added gamma_var
	
    
    size_t r, pos;
	for (size_t i=0; i<ns_test; ++i) {
        
        // save the data along the order of all variants, snp_pos is sorted by order
		outfile<<snp_pos[i].pos << "\t" << snp_pos[i].rs<<"\t" << snp_pos[i].chr<<"\t"
		<<snp_pos[i].bp<<"\t";
		
        r = mapOrder2Rank[i];
        pos = mapOrder2pos[i];
        if (pos == snp_pos[i].pos) {
      
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, i)<<"\t" << pos_loglr[r].second << "\t";
        
		if (beta_g[pos].second!=0) {
			outfile<<beta_g[pos].first/beta_g[pos].second<<"\t"<< beta_g[pos].second/(double)w <<endl;
		}
		else {
			outfile<<0.0<<"\t"<<0.0 <<endl;
		}
            
        }
        else {cerr << "pos dose not match snp_pos[i].pos...\n"; exit(1);}
        
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteParam (const gsl_vector *alpha) 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		

		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
				<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";				
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		outfile<<0.0<<"\t"<<0.0<<endl;
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const gsl_matrix *Result_theta, const gsl_matrix *Result_sigma, const size_t w_col)
{
	string file_gamma, file_hyp, file_theta, file_sigma;
	file_gamma="./output/"+file_out;
	file_gamma+=".gamma.txt";
	file_hyp="./output/"+file_out;
	file_hyp+=".hyp.txt";
    file_theta="./output/"+file_out;
	file_theta+=".theta.txt";
    file_sigma="./output/"+file_out;
	file_sigma+=".sigma.txt";

	ofstream outfile_gamma, outfile_hyp, outfile_theta, outfile_sigma;
		
	if (flag==0) {
		outfile_gamma.open (file_gamma.c_str(), ofstream::out);
		outfile_hyp.open (file_hyp.c_str(), ofstream::out);
        outfile_theta.open (file_theta.c_str(), ofstream::out);
		outfile_sigma.open (file_sigma.c_str(), ofstream::out);
        
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
        if (!outfile_theta) {cout<<"error writing file: "<<file_theta<<endl; return;}
		if (!outfile_sigma) {cout<<"error writing file: "<<file_sigma<<endl; return;}
        
		outfile_hyp<<"pve \t pge \t n_gamma \t logPosterior"<<endl;
		
		for (size_t i=0; i<s_max; ++i) {
			outfile_gamma<<"s"<<i<<"\t";
		}
		outfile_gamma << endl;
        
        for (size_t i=0; i < n_type; i++) {
            outfile_theta << "func_type"<<i << "\t";
            outfile_sigma << "func_type"<<i << "\t";
        }
        outfile_theta << endl;
        outfile_sigma << endl;
	}
	else {
		outfile_gamma.open (file_gamma.c_str(), ofstream::app);
		outfile_hyp.open (file_hyp.c_str(), ofstream::app);
        outfile_theta.open (file_theta.c_str(), ofstream::out);
		outfile_sigma.open (file_sigma.c_str(), ofstream::out);
        
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
        if (!outfile_theta) {cout<<"error writing file: "<<file_theta<<endl; return;}
		if (!outfile_sigma) {cout<<"error writing file: "<<file_sigma<<endl; return;}
        
		size_t w;
		if (w_col==0) {w=w_pace;}
		else {w=w_col;}
		
		for (size_t i=0; i<w; ++i) {
			outfile_hyp<<scientific;
			outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, 0)<<gsl_matrix_get (Result_hyp, i, 1)<<(int)gsl_matrix_get (Result_hyp, i, 2)<<"\t" << setprecision(6) << gsl_matrix_get (Result_hyp, i, 3)<<endl;
		}
		
		for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<s_max; ++j) {
				outfile_gamma<<(int)gsl_matrix_get (Result_gamma, i, j)<<"\t";
			}
			outfile_gamma<<endl;
		}
        
        for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<n_type; ++j) {
				outfile_theta<<scientific<<setprecision(6) << gsl_matrix_get (Result_theta, i, j)<<"\t";
			}
			outfile_theta<<endl;
		}
        
        for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<n_type; ++j) {
				outfile_sigma<<scientific<<setprecision(6)<< gsl_matrix_get (Result_sigma, i, j)<<"\t";
			}
			outfile_sigma<<endl;
		}
	}
	
	outfile_hyp.close();
	outfile_hyp.clear();
	outfile_gamma.close();
	outfile_gamma.clear();
    outfile_theta.close();
	outfile_theta.clear();
    outfile_sigma.close();
	outfile_sigma.clear();
	return;
}


//Revise p_gamma 
void BSLMM::CalcPgamma (double *p_gamma)
{
	double p, q;
    p = 0.9 / 1000;
    q = 0.1 / (ns_test-1000);
    
	for (size_t i=0; i<ns_test; ++i) {
        if(i < 1000) p_gamma[i] = p;
        else p_gamma[i] = q;
	}
	return;
}


void BSLMM::SetXgamma (gsl_matrix *Xgamma, uchar **X, vector<size_t> &rank)
{
	size_t pos;
	for (size_t i=0; i<rank.size(); ++i) {
		pos=mapRank2pos[rank[i]];        
		gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
        getGTgslVec(X, &Xgamma_col.vector, pos, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
    }	
	return;
}


void BSLMM::SetXgamma (uchar **X, const gsl_matrix *X_old, const gsl_matrix *XtX_old, const gsl_vector *Xty_old, const gsl_vector *y, const vector<size_t> &rank_old, const vector<size_t> &rank_new, gsl_matrix *X_new, gsl_matrix *XtX_new, gsl_vector *Xty_new)
{
    double d;
    // cout << "X_add set start" << endl;
    
    //rank_old and rank_new are sorted already inside PorposeGamma
    //calculate vectors rank_remove and rank_add
    //  size_t v_size=max(rank_old.size(), rank_new.size());
    //make sure that v_size is larger than repeat
    size_t v_size=20;
    vector<size_t> rank_remove(v_size), rank_add(v_size), rank_union(s_max+v_size);
    vector<size_t>::iterator it;
    
    it=set_difference (rank_old.begin(), rank_old.end(), rank_new.begin(), rank_new.end(), rank_remove.begin());
    rank_remove.resize(it-rank_remove.begin());
    
    it=set_difference (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_add.begin());
    rank_add.resize(it-rank_add.begin());
    
    it=set_union (rank_new.begin(), rank_new.end(), rank_old.begin(), rank_old.end(), rank_union.begin());
    rank_union.resize(it-rank_union.begin());
    
    //map rank_remove and rank_add
    map<size_t, int> mapRank2in_remove, mapRank2in_add;
    for (size_t i=0; i<rank_remove.size(); i++) {
        mapRank2in_remove[rank_remove[i]]=1;
    }
    for (size_t i=0; i<rank_add.size(); i++) {
        mapRank2in_add[rank_add[i]]=1;
    }
    
    //obtain the subset of matrix/vector
    gsl_matrix_const_view Xold_sub=gsl_matrix_const_submatrix(X_old, 0, 0, X_old->size1, rank_old.size());
    gsl_matrix_const_view XtXold_sub=gsl_matrix_const_submatrix(XtX_old, 0, 0, rank_old.size(), rank_old.size());
    gsl_vector_const_view Xtyold_sub=gsl_vector_const_subvector(Xty_old, 0, rank_old.size());
    
    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(X_new, 0, 0, X_new->size1, rank_new.size());
    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
    gsl_vector_view Xtynew_sub=gsl_vector_subvector(Xty_new, 0, rank_new.size());
    
    if (rank_remove.size()==0 && rank_add.size()==0) {
        gsl_matrix_memcpy(&Xnew_sub.matrix, &Xold_sub.matrix);
        gsl_matrix_memcpy(&XtXnew_sub.matrix, &XtXold_sub.matrix);
        gsl_vector_memcpy(&Xtynew_sub.vector, &Xtyold_sub.vector);
        //cout << "rank_old = rank_new; " << "Xgamma_new set success" << endl;
    } else {
        size_t i_old, j_old, i_new, j_new, i_add, j_add, i_flag, j_flag;
        if (rank_add.size()==0) {
            i_old=0; i_new=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}
                
                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                
                d=gsl_vector_get (Xty_old, i_old);
                gsl_vector_set (Xty_new, i_new, d);
                
                j_old=i_old; j_new=i_new;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}
                    
                    d=gsl_matrix_get(XtX_old, i_old, j_old);
                    
                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}
                    
                    j_old++; j_new++;
                }
                i_old++; i_new++;
            }
            //cout << "X_add = NULL; " << "Xgamma_new set success" << endl;
        } else {
            //rank_add has length > 0
            gsl_matrix *X_add=gsl_matrix_alloc(X_old->size1, rank_add.size() );
            gsl_matrix *XtX_aa=gsl_matrix_alloc(X_add->size2, X_add->size2);
            gsl_matrix *XtX_ao=gsl_matrix_alloc(X_add->size2, X_old->size2);
            gsl_vector *Xty_add=gsl_vector_alloc(X_add->size2);
            
            //get X_add
            SetXgamma (X_add, X, rank_add);
            //cout << "X_add set success" << endl;
            
            //get t(X_add)X_add and t(X_add)X_temp
            clock_t time_start=clock();
            
            //somehow the lapack_dgemm does not work here
            //#ifdef WITH_LAPACK
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_add, 0.0, XtX_aa);
            //lapack_dgemm ((char *)"T", (char *)"N", 1.0, X_add, X_old, 0.0, XtX_ao);
            
            //#else
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, X_add, 0.0, XtX_aa);
            gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, X_old, 0.0, XtX_ao);
            //#endif
            gsl_blas_dgemv(CblasTrans, 1.0, X_add, y, 0.0, Xty_add);
            
            time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            
            //save to X_new, XtX_new and Xty_new
            i_old=0; i_new=0; i_add=0;
            for (size_t i=0; i<rank_union.size(); i++) {
                if (mapRank2in_remove.count(rank_old[i_old])!=0) {i_old++; continue;}
                if (mapRank2in_add.count(rank_new[i_new])!=0) {i_flag=1; //within x_add
                } else {i_flag=0; //within x_common
                }
                
                gsl_vector_view Xnew_col=gsl_matrix_column(X_new, i_new);
                if (i_flag==1) {
                    gsl_vector_view Xcopy_col=gsl_matrix_column(X_add, i_add);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                } else {
                    gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(X_old, i_old);
                    gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
                }
                //  cout << "Xgamma_new set success" << endl;
                
                if (i_flag==1) {
                    d=gsl_vector_get (Xty_add, i_add);
                } else {
                    d=gsl_vector_get (Xty_old, i_old);
                }
                gsl_vector_set (Xty_new, i_new, d);
                // cout << "Xty_new set success" << endl;
                
                j_old=i_old; j_new=i_new; j_add=i_add;
                for (size_t j=i; j<rank_union.size(); j++) {
                    if (mapRank2in_remove.count(rank_old[j_old])!=0) {j_old++; continue;}
                    if (mapRank2in_add.count(rank_new[j_new])!=0) {j_flag=1;} else {j_flag=0;}
                    
                    if (i_flag==1 && j_flag==1) {
                        d=gsl_matrix_get(XtX_aa, i_add, j_add);          
                    } else if (i_flag==1) {
                        d=gsl_matrix_get(XtX_ao, i_add, j_old);
                    } else if (j_flag==1) {
                        d=gsl_matrix_get(XtX_ao, j_add, i_old);
                    } else {
                        d=gsl_matrix_get(XtX_old, i_old, j_old);
                    }
                    
                    gsl_matrix_set (XtX_new, i_new, j_new, d);
                    if (i_new!=j_new) {gsl_matrix_set (XtX_new, j_new, i_new, d);}
                    
                    j_new++; if (j_flag==1) {j_add++;} else {j_old++;}
                }
                //cout << "XtX_new success" << endl;
                i_new++; if (i_flag==1) {i_add++;} else {i_old++;}
            }
            // cout << "X_gamma set success" << endl;
            
            gsl_matrix_free(X_add);
            gsl_matrix_free(XtX_aa);
            gsl_matrix_free(XtX_ao);
            gsl_vector_free(Xty_add);
        }
        
    }
    
    rank_remove.clear();
    rank_add.clear();
    rank_union.clear();
    mapRank2in_remove.clear();
    mapRank2in_add.clear();
    
    return;
}


double BSLMM::CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2) 
{
	double pve, var_y;	
	
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *Xty=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *OiXty=gsl_vector_alloc (UtXgamma->size2);

	gsl_matrix_set_identity (Omega);
	gsl_matrix_scale (Omega, 1.0/sigma_a2); 

#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", 1.0, UtXgamma, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, UtXgamma, UtXgamma, 1.0, Omega);	
#endif
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma, Uty, 0.0, Xty);

	CholeskySolve(Omega, Xty, OiXty);
	
	gsl_blas_ddot (Xty, OiXty, &pve);
	gsl_blas_ddot (Uty, Uty, &var_y);
	
	pve/=var_y;
	
	gsl_matrix_free (Omega);
	gsl_vector_free (Xty);
	gsl_vector_free (OiXty);

	return pve;
}

bool comp_vec (size_t a, size_t b)
{
	return (a < b);
}

bool comp_lr (pair<size_t, double> a, pair<size_t, double> b)
{
	return (a.second > b.second);
}

//JY
bool comp_snp(const SNPPOS& lhs, const SNPPOS& rhs){
    return (lhs.chr.compare(rhs.chr) < 0) || ((lhs.chr.compare(rhs.chr) == 0) && (lhs.bp < rhs.bp));
}
//JY


void BSLMM::InitialMCMC (uchar **X, const gsl_vector *Uty, LModel &model, vector<pair<size_t, double> > &pos_loglr, const vector<snpPos> &snp_pos)

 {
 vector<pair<size_t, double> > rank_loglr;
 size_t posr, radd;
 
 double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
 model.cHyp.n_gamma=0;
 for (size_t i=0; i<pos_loglr.size(); ++i) {
     if (2.0*pos_loglr[i].second>q_genome) {model.cHyp.n_gamma++;}
 }
 if (model.cHyp.n_gamma<10) {model.cHyp.n_gamma=10;}
 if (model.cHyp.n_gamma>s_max) {model.cHyp.n_gamma=s_max;}
 if (model.cHyp.n_gamma<s_min) {model.cHyp.n_gamma=s_min;}
 cout << "number of snps = " << model.cHyp.n_gamma << endl;
 
 for (size_t i=1; i<1000; ++i) {
     rank_loglr.push_back(make_pair(i, pos_loglr[i].second));
 }
 cout << endl;
 
 model.rank.clear();
 model.rank.push_back(0);
 posr = mapRank2pos[0];
 
 gsl_matrix * Xr = gsl_matrix_alloc(ni_test, model.cHyp.n_gamma);
 gsl_matrix_set_zero(Xr);
 gsl_vector * xvec = gsl_vector_alloc(ni_test);
     getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
 
 gsl_matrix * XtXr = gsl_matrix_alloc(model.cHyp.n_gamma, model.cHyp.n_gamma);
 gsl_matrix_set_zero(XtXr);
 gsl_vector * Xtyr = gsl_vector_alloc(model.cHyp.n_gamma);
 gsl_vector_set_zero(Xtyr);
 
 gsl_vector * yres = gsl_vector_alloc(ni_test);
 gsl_vector * Xtxvec = gsl_vector_alloc(model.cHyp.n_gamma);
     gsl_vector_set_zero(Xtxvec);

 double xty, yty;
 gsl_blas_ddot(Uty, Uty, &yty);
 
 for (size_t i=1; i < model.cHyp.n_gamma; ++i){
     //cout << "i = " << i << "," ;
 gsl_matrix_set_col(Xr, (i-1), xvec);
 gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
 
 gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
 gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);
 
 gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, (i-1), 0, i);
 gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, (i-1), 0, i);
 gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
 gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);
 
 gsl_blas_ddot(xvec, Uty, &xty);
 gsl_vector_set(Xtyr, (i-1), xty);
 
 CalcRes(Xr, Uty, XtXr, Xtyr, yres, i, yty);
 for (size_t j=0; j<rank_loglr.size(); ++j) {
     posr = mapRank2pos[rank_loglr[j].first];
     getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
     rank_loglr[j].second = CalcLR(yres, xvec);
 }
 stable_sort (rank_loglr.begin(), rank_loglr.end(), comp_lr); 
 
 radd = rank_loglr[0].first;
 posr = mapRank2pos[radd];
     getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
 model.rank.push_back(radd);
 rank_loglr.erase(rank_loglr.begin());
 }
     stable_sort (model.rank.begin(), model.rank.end(), comp_vec); //sort the initial rank.
     PrintVector(model.rank);
     
 gsl_vector_free(Xtxvec);
 gsl_matrix_free(Xr);
 gsl_matrix_free(XtXr);
 gsl_vector_free(Xtyr);
 gsl_vector_free(xvec);
 gsl_vector_free(yres);
     
     
 vector<string> iniRank; //JY added vector iniRank to save all significant snp ID
 size_t order;
 for (size_t i=0; i<model.rank.size(); ++i) {
 order = mapRank2Order[model.rank[i]];
     
 iniRank.push_back(snp_pos[order].rs);
 }
 WriteIniRank(iniRank);  // write out initial sig snp ID
 
 model.cHyp.logp=log((double)model.cHyp.n_gamma/(double)ns_test);
 model.cHyp.h=pve_null;
 
 if (model.cHyp.logp==0) {model.cHyp.logp=-0.000001;}
 if (model.cHyp.h==0) {model.cHyp.h=0.1;}
 
 gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, model.cHyp.n_gamma);
 SetXgamma (UtXgamma, X, model.rank);
 double sigma_a2;
 if (trace_G!=0) {
 sigma_a2=model.cHyp.h/(trace_G*(1-model.cHyp.h)*exp(model.cHyp.logp)*(double)ns_test);
     cout << "Initial sigma_a2 = " << sigma_a2 << endl;
 } else {
 sigma_a2=model.cHyp.h/( (1-model.cHyp.h)*exp(model.cHyp.logp)*(double)ns_test);
 }
 if (sigma_a2==0) {sigma_a2=0.025;}
 model.cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/model.cHyp.h;
 gsl_matrix_free (UtXgamma);
 
 if (model.cHyp.rho>1.0) {model.cHyp.rho=1.0;}
 
 if (model.cHyp.h<h_min) {model.cHyp.h=h_min;}
 if (model.cHyp.h>h_max) {model.cHyp.h=h_max;}
 if (model.cHyp.rho<rho_min) {model.cHyp.rho=rho_min;}
 if (model.cHyp.rho>rho_max) {model.cHyp.rho=rho_max;}
 if (model.cHyp.logp<logp_min) {model.cHyp.logp=logp_min;}
 if (model.cHyp.logp>logp_max) {model.cHyp.logp=logp_max;}
 
 cout<<"initial value of h = "<<model.cHyp.h<<endl;
 cout<<"initial value of rho = "<<model.cHyp.rho<<endl;
 cout<<"initial value of pi = "<<exp(model.cHyp.logp)<<endl;
 cout<<"initial value of |gamma| = "<<model.cHyp.n_gamma<<endl;
 
 return;
 }

void BSLMM::InitialMCMC (uchar **X, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<snpPos> &snp_pos)

{
    vector<pair<size_t, double> > rank_loglr;
    size_t posr, radd;
    
    double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
    cHyp.n_gamma=0;
    for (size_t i=0; i<pos_loglr.size(); ++i) {
        if (2.0*pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
    }
    if (cHyp.n_gamma<10) {cHyp.n_gamma=10;}
    if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
    if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}
    cout << "number of snps = " << cHyp.n_gamma << endl;
    
    for (size_t i=1; i<1000; ++i) {
        rank_loglr.push_back(make_pair(i, pos_loglr[i].second));
    }
    cout << endl;
    
    rank.clear();
    rank.push_back(0);
    posr = mapRank2pos[0];
   // cout << "rank added: " << 0 << ", ";

    
    gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
    gsl_matrix_set_zero(Xr);
    gsl_vector * xvec = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize); //get geno column
    
    gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
    gsl_matrix_set_zero(XtXr);
    gsl_vector * Xtyr = gsl_vector_alloc(cHyp.n_gamma);
    gsl_vector_set_zero(Xtyr);
    
    gsl_vector * yres = gsl_vector_alloc(ni_test);
    gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);
    gsl_vector_set_zero(Xtxvec);
    
    double xty, yty;
    gsl_blas_ddot(Uty, Uty, &yty);
    
    for (size_t i=1; i < cHyp.n_gamma; ++i){
        gsl_matrix_set_col(Xr, (i-1), xvec);
        gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
        
        gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
        gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);
        
        gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, (i-1), 0, i);
        gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, (i-1), 0, i);
        gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
        gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);
        
        gsl_blas_ddot(xvec, Uty, &xty);
        gsl_vector_set(Xtyr, (i-1), xty);
        
        CalcRes(Xr, Uty, XtXr, Xtyr, yres, i, yty);
        for (size_t j=0; j<rank_loglr.size(); ++j) {
            posr = mapRank2pos[rank_loglr[j].first];
            getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize); // get geno column
            rank_loglr[j].second = CalcLR(yres, xvec);
        }
        stable_sort (rank_loglr.begin(), rank_loglr.end(), comp_lr); //sort the initial rank.
        
        radd = rank_loglr[0].first;
       // cout << "rank added: " << radd << ", ";
        posr = mapRank2pos[radd];
        getGTgslVec(X, xvec, posr, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
        rank.push_back(radd);
        rank_loglr.erase(rank_loglr.begin());
    }
    
    gsl_matrix_free(Xr);
    gsl_matrix_free(XtXr);
    gsl_vector_free(Xtyr);
    gsl_vector_free(xvec);
    gsl_vector_free(yres);
    gsl_vector_free(Xtxvec);

    stable_sort (rank.begin(), rank.end(), comp_vec); //sort the initial rank.
    PrintVector(rank);
    
    vector<string> iniRank; //JY added vector iniRank to save all significant snp ID
    size_t order;
    for (size_t i=0; i<rank.size(); ++i) {
        order = mapRank2Order[rank[i]];
        iniRank.push_back(snp_pos[order].rs);
    }
    WriteIniRank(iniRank);  // write out initial sig snp ID
    
    cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
    cHyp.h=pve_null;
    
    if (cHyp.logp==0) {cHyp.logp=-0.000001;}
    if (cHyp.h==0) {cHyp.h=0.1;}
    
    gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
    SetXgamma (UtXgamma, X, rank);
    double sigma_a2;
    if (trace_G!=0) {
        sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
    } else {
        sigma_a2=cHyp.h*1.0/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
    }
    if (sigma_a2==0) {sigma_a2=0.025;}
    cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
    gsl_matrix_free (UtXgamma);
    
    if (cHyp.rho>1.0) {cHyp.rho=1.0;}
    
    if (cHyp.h<h_min) {cHyp.h=h_min;}
    if (cHyp.h>h_max) {cHyp.h=h_max;}
    if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
    if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
    if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
    if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}
    
    cout<<"initial value of h = "<<cHyp.h<<endl;
    cout<<"initial value of rho = "<<cHyp.rho<<endl;
    cout<<"initial value of pi = "<<exp(cHyp.logp)<<endl;
    cout<<"initial value of |gamma| = "<<cHyp.n_gamma<<endl;
    
    return;
}


/*
void BSLMM::InitialMCMC (const gsl_matrix *UtX, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<snpPos> &snp_pos)
{
	double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
	
	cHyp.n_gamma=0;
	for (size_t i=0; i<pos_loglr.size(); ++i) {
		if (2.0*pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
	}
	if (cHyp.n_gamma<10) {cHyp.n_gamma=10;}
	
	if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
	if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}	
	
    vector<string> iniRank; //JY added vector iniRank to save all significant snp ID
    size_t order;
    
	rank.clear();
	for (size_t i=0; i<cHyp.n_gamma; ++i) {
		rank.push_back(i);
        order = mapRank2Order[i];
        iniRank.push_back(snp_pos[order].rs);
	}
    //WriteIniRank(iniRank); // write out initial sig snp ID
    
	
	cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
	cHyp.h=pve_null; 
	
	if (cHyp.logp==0) {cHyp.logp=-0.000001;}
	if (cHyp.h==0) {cHyp.h=0.1;}

	gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
	SetXgamma (UtXgamma, UtX, rank);
	double sigma_a2;
	if (trace_G!=0) {
	  sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	} else {
	  sigma_a2=cHyp.h*1.0/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	}
	if (sigma_a2==0) {sigma_a2=0.025;}	
	cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
	gsl_matrix_free (UtXgamma);
	
	if (cHyp.rho>1.0) {cHyp.rho=1.0;}
	
	if (cHyp.h<h_min) {cHyp.h=h_min;}
	if (cHyp.h>h_max) {cHyp.h=h_max;}
	if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
	if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
	if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
	if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}
	
	
//	if (fix_sigma>=0) {
//		fix_sigma=cHyp.h;
//		rho_max=1-cHyp.h;
//		cHyp.rho=rho_max/2.0;
//	}
	
	//Initial for grid sampling:
//	cHyp.h=0.225;
//	cHyp.rho=1.0;
//	cHyp.logp=-4.835429;
	
	cout<<"initial value of h = "<<cHyp.h<<endl;
	cout<<"initial value of rho = "<<cHyp.rho<<endl;
	cout<<"initial value of pi = "<<exp(cHyp.logp)<<endl;
	cout<<"initial value of |gamma| = "<<cHyp.n_gamma<<endl;
	
	return;
}*/


double BSLMM::CalcPosterior (const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *Utu, gsl_vector *alpha_prime, class HYPBSLMM &cHyp)
{
	double sigma_b2=cHyp.h*(1.0-cHyp.rho)/(trace_G*(1-cHyp.h));
	
	gsl_vector *Utu_rand=gsl_vector_alloc (Uty->size);	
	gsl_vector *weight_Hi=gsl_vector_alloc (Uty->size);
	
	double logpost=0.0;
	double d, ds, uy, Hi_yy=0, logdet_H=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		ds=d/(d+1.0);
		d=1.0/(d+1.0);		
		gsl_vector_set (weight_Hi, i, d);
		
		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		Hi_yy+=d*uy*uy;
		
		gsl_vector_set (Utu_rand, i, gsl_ran_gaussian(gsl_r, 1)*sqrt(ds));
	}
	
	//sample tau
	double tau=1.0;
	if (a_mode==11) {tau = gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/Hi_yy); }
	
	//sample alpha
	gsl_vector_memcpy (alpha_prime, Uty);
	gsl_vector_mul (alpha_prime, weight_Hi);
	gsl_vector_scale (alpha_prime, sigma_b2);
	
	//sample u
	gsl_vector_memcpy (Utu, alpha_prime);
	gsl_vector_mul (Utu, K_eval);
	if (a_mode==11) {gsl_vector_scale (Utu_rand, sqrt(1.0/tau));}
	gsl_vector_add (Utu, Utu_rand);	
	
	//for quantitative traits, calculate pve and ppe
	if (a_mode==11) {
		gsl_blas_ddot (Utu, Utu, &d);
		cHyp.pve=d/(double)ni_test;	
		cHyp.pve/=cHyp.pve+1.0/tau;
		cHyp.pge=0.0;	
	}

	//calculate likelihood
	logpost=-0.5*logdet_H;
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(Hi_yy);}
	else {logpost-=0.5*Hi_yy;}
	
	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1-exp(cHyp.logp));
	
	gsl_vector_free (Utu_rand);
	gsl_vector_free (weight_Hi);
	
	return logpost;
}


double BSLMM::CalcPosterior (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *UtXb, gsl_vector *Utu, gsl_vector *alpha_prime, gsl_vector *beta, class HYPBSLMM &cHyp)
{
	clock_t time_start;	
	
	double sigma_a2=cHyp.h*cHyp.rho/(trace_G*(1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	double sigma_b2=cHyp.h*(1.0-cHyp.rho)/(trace_G*(1-cHyp.h));
	
	double logpost=0.0;
	double d, ds, uy, P_yy=0, logdet_O=0.0, logdet_H=0.0;
	
	gsl_matrix *UtXgamma_eval=gsl_matrix_alloc (UtXgamma->size1, UtXgamma->size2);	
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *XtHiy=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *beta_hat=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *Utu_rand=gsl_vector_alloc (UtXgamma->size1);	
	gsl_vector *weight_Hi=gsl_vector_alloc (UtXgamma->size1);
	
	gsl_matrix_memcpy (UtXgamma_eval, UtXgamma);
	
	logdet_H=0.0; P_yy=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		gsl_vector_view UtXgamma_row=gsl_matrix_row (UtXgamma_eval, i);
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		ds=d/(d+1.0);
		d=1.0/(d+1.0);
		gsl_vector_set (weight_Hi, i, d);
		
		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		P_yy+=d*uy*uy;
		gsl_vector_scale (&UtXgamma_row.vector, d);
		
		gsl_vector_set (Utu_rand, i, gsl_ran_gaussian(gsl_r, 1)*sqrt(ds));
	}
	
	//calculate Omega
	gsl_matrix_set_identity (Omega);
	
	time_start=clock();
#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#endif	
	time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
	
	
	//calculate beta_hat
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma_eval, Uty, 0.0, XtHiy);	

	logdet_O=CholeskySolve(Omega, XtHiy, beta_hat);
	
	gsl_vector_scale (beta_hat, sigma_a2);

	gsl_blas_ddot (XtHiy, beta_hat, &d);
	P_yy-=d;
	
	//sample tau
	double tau=1.0;
	if (a_mode==11) {tau =gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/P_yy); }

	//sample beta
	for (size_t i=0; i<beta->size; i++)
	{
		d=gsl_ran_gaussian(gsl_r, 1); 
		gsl_vector_set(beta, i, d); 
	}
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, beta); 
	
	
	//it compuates inv(L^T(Omega)) %*% beta;  
	gsl_vector_scale(beta, sqrt(sigma_a2/tau));
	gsl_vector_add(beta, beta_hat); 
	gsl_blas_dgemv (CblasNoTrans, 1.0, UtXgamma, beta, 0.0, UtXb);
	
	//sample alpha
	gsl_vector_memcpy (alpha_prime, Uty);
	gsl_vector_sub (alpha_prime, UtXb);
	gsl_vector_mul (alpha_prime, weight_Hi);
	gsl_vector_scale (alpha_prime, sigma_b2);
	
	//sample u
	gsl_vector_memcpy (Utu, alpha_prime);
	gsl_vector_mul (Utu, K_eval);
	
	if (a_mode==11) {gsl_vector_scale (Utu_rand, sqrt(1.0/tau));}
	gsl_vector_add (Utu, Utu_rand);	
	
	
	//for quantitative traits, calculate pve and pge
	if (a_mode==11) {
		gsl_blas_ddot (UtXb, UtXb, &d);
		cHyp.pge=d/(double)ni_test;
	
		gsl_blas_ddot (Utu, Utu, &d);
		cHyp.pve=cHyp.pge+d/(double)ni_test;
		
		if (cHyp.pve==0) {cHyp.pge=0.0;}
		else {cHyp.pge/=cHyp.pve;}
		cHyp.pve/=cHyp.pve+1.0/tau;	
	}	
	

	gsl_matrix_free (UtXgamma_eval);
	gsl_matrix_free (Omega);
	gsl_vector_free (XtHiy);
	gsl_vector_free (beta_hat);
	gsl_vector_free (Utu_rand);	
	gsl_vector_free (weight_Hi);
	
	logpost=-0.5*logdet_H-0.5*logdet_O;
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(P_yy);}
	else {logpost-=0.5*P_yy;}
//	else {logpost+=-0.5*P_yy*tau+0.5*(double)ni_test*log(tau);}
	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1.0-exp(cHyp.logp));
	
	return logpost;
}

// from previous version
double BSLMM::CalcPosterior (const double yty, class HYPBSLMM &cHyp, const gsl_vector *pi_vec, const gsl_vector *sigma_vec, const vector<size_t> &rank)
{
    double logpost=0.0;
    
    //for quantitative traits, calculate pve and pge
    //pve and pge for case/control data are calculted in CalcCC_PVEnZ
    if (a_mode==11) {
        cHyp.pve=0.0;
        cHyp.pge=1.0;
    }
    
    //calculate likelihood
    if (a_mode==11) {logpost-=0.5*(double)ni_test*log(yty);}
    else {logpost-=0.5*yty;}
    
    //Calc prior portion
    double pi_rank=0.0, pi_nonrank=0.0, sumlogV = 0.0, sumV = 0.0, sumTheta = 0.0;
    for (size_t i=0; i < ns_test; i++) {
        pi_nonrank += log(1.0 - gsl_vector_get(pi_vec, i));
    }
    for (size_t i=0; i < s_size; i++) {
        order_i = mapRank2Order[rank[i]];
        pi_rank += log(gsl_vector_get(pi_vec, order_i));
        pi_nonrank -= log(1.0 - gsl_vector_get(pi_vec, order_i));
        
        sumV += 1.0 / gsl_vector_get(sigma_vec, i);
        sumlogV += log(gsl_vector_get(sigma_vec, i));
    }
    for (size_t i=0; i < n_type; i++) {
        sumTheta -= cHyp.log_theta[i];
    }
    cout << "pi_rank = " << pi_rank << "; pi_nonrank = " << pi_nonrank << endl;
    cout << "sumV = " << sumV << "; sumlogV = " << sumlogV << endl;
    
    // calculate posterior
    logpost += (1.0 - e) * sumlogV - e * sumV + pi_rank + pi_nonrank + sumTheta;
    
    //logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1-exp(cHyp.logp));
    
    return logpost;
}

void BSLMM::CalcVPi(class HYPBSLMM &cHyp, const gsl_vector *pi_vec, const gsl_vector *sigma_vec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos){
    
    if (cHyp.n_gamma != rank.size()) {
        cerr << "Error: cHyp.n_gamma not equal to the size of rank"
    }
    
    size_t order_i;
    vector<double> theta, weight_temp;
    for (size_t j=0; j < n_type; j++) {
        theta.push_back(exp(cHyp.log_theta[j]));
    }
    gsl_vector_set_zero(pi_vec);
    gsl_vector_set_zero(sigma_vec);
    double pi_temp, sigma_temp;
    
    for (size_t i=0; i < ns_test; i++) {
        weight_temp = snp_pos[order_i].weight;
        pi_temp = 0.0;
        for (size_t j=0; j < n_type; j++) {
            pi_temp += theta[j] * weight_temp[j];
        }
        gsl_vector_set(pi_vec, i, pi_temp);
    }
    
    for (size_t i=0; i < cHyp.n_gamma; i++) {
        order_i = mapRank2Order[rank[i]];
        weight_temp = snp_pos[order_i].weight;
        sigma_temp = 0.0;        
        for (size_t j=0; j < n_type; j++) {
            sigma_temp += cHyp.subvar[j] * weight_temp[j];
        }
        gsl_vector_set(sigma_vec, i, sigma_temp);
    }

}

double BSLMM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *pi_vec, gsl_vector *sigma_vec, const vector<size_t> &rank)
{
    //double sigma_a2=cHyp.h/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test * trace_G);
    
    double logpost=0.0;
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = rank.size(), order_i;
    
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *Omega_temp=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *V_beta=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);
    gsl_vector *beta_hat_temp=gsl_vector_alloc (s_size);
    gsl_vector *Xty_temp=gsl_vector_alloc (s_size);
        
    //calculate Omega
    gsl_matrix_set_identity (V_beta);
    gsl_vector_view V_beta_diag = gsl_matrix_diagonal(V_beta);
    gsl_vector_mul(&V_beta_diag.vector, &sigma_sub.vector); //inv(V_beta)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &XtX_sub.matrix, V_beta, 0.0, Omega);
    gsl_matrix_set_identity (M_temp);
    gsl_matrix_add (Omega, M_temp);
    
    //calculate beta_hat
    gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
    logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat_temp);	//solve Omega * beta_hat = Xty for beta_hat
    // Omega was inverted here
    // logdet_0 = det(Omega)
    gsl_blas_dgemv(CblasNoTrans, 1.0, V_beta, beta_hat_temp, 0.0, beta_hat);
    gsl_blas_ddot (Xty_temp, beta_hat, &d);
    P_yy-=d;
    
    //sample tau
    double tau=1.0;
    if (a_mode==11) {tau =gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/P_yy); }
    
    //sample beta
    for (size_t i=0; i<s_size; i++)
    {
        d=gsl_ran_gaussian(gsl_r, 1);
        gsl_vector_set(beta, i, d);
    }
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Omega, V_beta, 0.0, Omega_temp);
    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega_temp, &beta_sub.vector);
    
    gsl_vector_scale(&beta_sub.vector, sqrt(tau));
    gsl_vector_add(&beta_sub.vector, beta_hat);
    gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);
    
    //for quantitative traits, calculate pve and pge
    if (a_mode==11) {
        gsl_blas_ddot (Xb, Xb, &d);
        cHyp.pve=d/(double)ni_test;
        cHyp.pve/=cHyp.pve+1.0/tau;
        cHyp.pge=1.0;
    }
    
    logpost=-0.5*logdet_O;
    if (a_mode==11) {logpost-=0.5*(double)ni_test*log(P_yy);}
    else {logpost-=0.5*P_yy;}
    
    //Calc prior portion for posterior
    double pi_rank=0.0, pi_nonrank=0.0, sumlogV = 0.0, sumV = 0.0, sumTheta = 0.0;
    for (size_t i=0; i < ns_test; i++) {
         pi_nonrank += log(1.0 - gsl_vector_get(pi_vec, i));
    }
    for (size_t i=0; i < s_size; i++) {
        order_i = mapRank2Order[rank[i]];
        pi_rank += log(gsl_vector_get(pi_vec, order_i));
        pi_nonrank -= log(1.0 - gsl_vector_get(pi_vec, order_i));
        
        sumV += 1.0 / gsl_vector_get(sigma_vec, i);
        sumlogV += log(gsl_vector_get(sigma_vec, i));
    }
    for (size_t i=0; i < n_type; i++) {
        sumTheta -= cHyp.log_theta[i];
    }
    cout << "pi_rank = " << pi_rank << "; pi_nonrank = " << pi_nonrank << endl;
    cout << "sumV = " << sumV << "; sumlogV = " << sumlogV << endl;

    logpost += (1.0 - e) * sumlogV - e * sumV + pi_rank + pi_nonrank + sumTheta;
    //logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1.0-exp(cHyp.logp));
    
    gsl_matrix_free (Omega);
    gsl_matrix_free (Omega_temp);
    gsl_matrix_free (V_beta);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (beta_hat_temp);
    gsl_vector_free (Xty_temp);

    return logpost;
}

double BSLMM::ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, uchar **X, const gsl_vector *z, const gsl_matrix *Xgamma_old, const gsl_matrix *XtX_old, const gsl_vector *Xtz_old, const double &ztz, int &flag_gamma)
{
    map<size_t, int> mapRank2in;
    double unif, logp = 0.0;
    //int flag_gamma;
    size_t r_add, r_remove, col_id, r;
    
    rank_new.clear();
    if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"size wrong"<<endl;}
    
    if (cHyp_old.n_gamma!=0) {
        for (size_t i=0; i<rank_old.size(); ++i) {
            r=rank_old[i];
            rank_new.push_back(r);
            mapRank2in[r]=1;
        }
    }
    cHyp_new.n_gamma=cHyp_old.n_gamma;
    
    for (size_t i=0; i<repeat; ++i) {
        
        unif=gsl_rng_uniform(gsl_r);
        
        if (unif < 0.33 && cHyp_new.n_gamma<s_max) {flag_gamma=1;}
        else if (unif>=0.33 && unif < 0.67 && cHyp_new.n_gamma>s_min) {flag_gamma=2;}
        else if (unif>=0.67 && cHyp_new.n_gamma>0 && cHyp_new.n_gamma<ns_test) {flag_gamma=3;}
        else {flag_gamma=4;}
        
        if(flag_gamma==1)  {//add a snp;
            //cout << "add a snp" << endl;
            do {
                r_add=gsl_ran_discrete (gsl_r, gsl_t);
            } while (mapRank2in.count(r_add)!=0);
            
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            
            mapRank2in[r_add]=1;
            rank_new.push_back(r_add);
            cHyp_new.n_gamma++;
            logp+=-log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);
            // cout << "succesfully added a snp" << endl;
            
        }
        else if (flag_gamma==2) {//delete a snp;
            // cout << "delete a snp" << endl;
            
            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];
            
            double prob_total=1.0;
            for (size_t ii=0; ii<cHyp_new.n_gamma; ++ii) {
                r=rank_new[ii];
                prob_total-=p_gamma[r];
            }
            prob_total+=p_gamma[r_remove];
            
            mapRank2in.erase(r_remove);
            rank_new.erase(rank_new.begin()+col_id);
            logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
            cHyp_new.n_gamma--;
            //  cout << "succesfully deleted a snp" << endl;
        }
        else if (flag_gamma==3) {//switch a snp;
            // cout << "switch a snp" << endl;
            long int o_add, o_remove;
            long int o_rj, o_aj;
            size_t j_add, j_remove, o;
            
            gsl_matrix *Xgamma_temp=gsl_matrix_alloc (ni_test, s_max);
            gsl_matrix *XtX_gamma=gsl_matrix_alloc (s_max, s_max);
            gsl_vector *Xtz_gamma=gsl_vector_alloc (s_max);
            gsl_vector *z_res = gsl_vector_alloc(ni_test);
            gsl_ran_discrete_t *gsl_s, *gsl_a; //JY added dynamic gsl_s
            
            double *p_BFr = new double[ns_neib];
            double *p_BFa = new double[ns_neib];
            
            col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
            r_remove=rank_new[col_id];//careful with the proposal
            if(mapRank2in.count(r_remove) == 0) {cout << "wrong proposal of r_remove;" << endl; exit(1);}
            o_remove = mapRank2Order[r_remove];
            rank_new.erase(rank_new.begin()+col_id);
            size_t s_size = rank_new.size();
            mapRank2in.erase(r_remove);
            
            if (cHyp_new.n_gamma<=20 || cHyp_old.n_gamma<=20) {
                SetXgamma (Xgamma_temp, X, rank_new);
                CalcXtX (Xgamma_temp, z, s_size, XtX_gamma, Xtz_gamma);
            } else {
                SetXgamma (X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_temp, XtX_gamma, Xtz_gamma);
            }
            
            CalcRes(Xgamma_temp, z, XtX_gamma, Xtz_gamma, z_res, s_size, ztz);
            gsl_s = MakeProposal(o_remove, p_BFr, X, z_res, mapRank2in);
            //construct gsl_s, JY
            
            j_add = gsl_ran_discrete(gsl_r, gsl_s);
            o_add = (o_remove - win) + j_add;
            r_add = mapOrder2Rank[o_add];
            //cout << "o_add = " << o_add <<  "; r_add = "<<r_add << endl;
            if((o_add < 0) || (o_add >= (long int)ns_test) || (o_add == (long int)o_remove))
                cout << "ERROR proposing switch snp"; //new snp != removed snp
            
            gsl_a = MakeProposal(o_add, p_BFa, X, z_res, mapRank2in);
            
            double prob_total_remove=1.0;
            double prob_total_add=1.0;
            
            for (size_t ii=0; ii<rank_new.size(); ++ii) {
                r = rank_new[ii];
                o = mapRank2Order[r];
                o_rj = ((long int)o - o_remove) + win;
                o_aj = ((long int)o - o_add) + win;
                if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
                if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
            }
            
            j_remove = o_remove - o_add + win;
            logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
            logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)
            
            mapRank2in[r_add]=1;
            rank_new.push_back(r_add);
            
            gsl_matrix_free(Xgamma_temp);
            gsl_matrix_free(XtX_gamma);
            gsl_vector_free(Xtz_gamma);
            gsl_vector_free(z_res);
            gsl_ran_discrete_free(gsl_s);
            gsl_ran_discrete_free(gsl_a);
            
            delete[] p_BFr;
            delete[] p_BFa;
            //cout << "successfully switched a snp" << endl;
        }
        
        else {logp+=0.0;}//do not change
        
    }
    
    stable_sort (rank_new.begin(), rank_new.end(), comp_vec);
    mapRank2in.clear();
    return logp;
}

void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col, const vector<snpPos> &snp_pos, const vector<pair<size_t, double> > &pos_loglr)
{
    string file_gamma, file_hyp;
    file_gamma="./output/"+file_out;
    file_gamma+=".gamma.txt";
    file_hyp="./output/"+file_out;
    file_hyp+=".hyp.txt";

    
    ofstream outfile_gamma, outfile_hyp;
    
    if (flag==0) {
        outfile_gamma.open (file_gamma.c_str(), ofstream::out);
        outfile_hyp.open (file_hyp.c_str(), ofstream::out);
        
        if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
        if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
        
        outfile_hyp<<"h \t pve \t rho \t pge \t pi \t n_gamma"<<endl;
        
        for (size_t i=0; i<s_max; ++i) {
            outfile_gamma<<"s"<<i<<"\t";
        }
        outfile_gamma << endl;
    }
    else {
        outfile_gamma.open (file_gamma.c_str(), ofstream::app);
        outfile_hyp.open (file_hyp.c_str(), ofstream::app);
        
        if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
        if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
        
        size_t w;
        if (w_col==0) {w=w_pace;}
        else {w=w_col;}
        
        for (size_t i=0; i<w; ++i) {
            outfile_hyp<<scientific;
            for (size_t j=0; j<4; ++j) {
                outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, j)<<"\t";
            }
            outfile_hyp<<setprecision(6)<<exp(gsl_matrix_get (Result_hyp, i, 4))<<"\t";
            outfile_hyp<<(int)gsl_matrix_get (Result_hyp, i, 5)<<"\t";
            outfile_hyp<<endl;
        }
        
        for (size_t i=0; i<w; ++i) {
            for (size_t j=0; j<s_max; ++j) {
                outfile_gamma<<(int)gsl_matrix_get (Result_gamma, i, j)<<"\t";
            }
            outfile_gamma<<endl;
        }
        
    }
    
    outfile_hyp.close();
    outfile_hyp.clear();
    outfile_gamma.close();
    outfile_gamma.clear();

    return;
}

void BSLMM::WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".param.txt";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    outfile<<"chr"<<"\t"<<"rs"<<"\t"
    <<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
    <<"beta"<<"\t"<<"gamma" << endl; //JY added gamma_var
    
    size_t t=0;
    for (size_t i=0; i<ns_total; ++i) {
        if (indicator_snp[i]==0) {continue;}
        
        outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
        <<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";
        
        outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
        if (beta_g[t].second!=0) {
            outfile<<beta_g[t].first/beta_g[t].second<<"\t"<< beta_g[t].second/(double)w <<endl;
        }
        else {
            outfile<<0.0<<"\t"<<0.0 <<endl;
        }
        t++;
    }		
    
    outfile.clear();	
    outfile.close();	
    return;
}


void BSLMM::MCMC (uchar **X, const gsl_vector *y, bool original_method) {
    
    if (original_method) {
        cout << "Run previous working version of MCMC...\n";
    }
    
    clock_t time_start;
    
    double time_set=0, time_post=0;
    
    class HYPBSLMM cHyp_old, cHyp_new;
    
    //new model related
    gsl_vector *pi_vec = gsl_vector_alloc (ns_test);
    gsl_vector *sigma_vec = gsl_vector_alloc(s_max);
    gsl_matrix *Result_theta = gsl_matrix_alloc(w_pace, n_type); // save theta samples
    gsl_matrix *Result_sigma = gsl_matrix_alloc(w_pace, n_type); // save sigma samples
    
    //same as old model
    gsl_matrix *Result_hyp=gsl_matrix_alloc (w_pace, 4);
    gsl_matrix *Result_gamma=gsl_matrix_alloc (w_pace, s_max);
    
    gsl_vector *Xb_new=gsl_vector_alloc (ni_test);
    gsl_vector *Xb_old=gsl_vector_alloc (ni_test);
    gsl_vector *z_hat=gsl_vector_alloc (ni_test);
    gsl_vector *z=gsl_vector_alloc (ni_test);
    
    gsl_matrix *Xgamma_old=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_old=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_old=gsl_vector_alloc (s_max);
    gsl_vector *beta_old=gsl_vector_alloc (s_max);
    
    gsl_matrix *Xgamma_new=gsl_matrix_alloc (ni_test, s_max);
    gsl_matrix *XtX_new=gsl_matrix_alloc (s_max, s_max);
    gsl_vector *Xtz_new=gsl_vector_alloc (s_max);
    gsl_vector *beta_new=gsl_vector_alloc (s_max);
    
    double ztz=0.0;
    gsl_vector_memcpy (z, y);
    //for quantitative traits, y is centered already in gemma.cpp, but just in case
    double mean_z = CenterVector (z);
    
    gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    cout << "ztz = " << ztz << endl;
    cout << "mean of z = " << mean_z << endl;
    
    double logPost_new, logPost_old;
    double logMHratio;
    
    gsl_matrix_set_zero (Result_gamma);
    if (a_mode==13) {
        pheno_mean=0.0;
    }
    
    vector<pair<double, double> > beta_g;
    for (size_t i=0; i<ns_test; i++) {
        beta_g.push_back(make_pair(0.0, 0.0));
    }
    
    vector<size_t> rank_new, rank_old;
    vector<pair<size_t, double> > pos_loglr;
    
    time_start=clock();
    MatrixCalcLmLR (X, z, pos_loglr, ns_test, ni_test, trace_G, CompBuffSizeVec, UnCompBufferSize); //const genotype
    cout << "trace_G = trace(X'X) = " << trace_G << endl;
    
    // Jingjing add a vector of "snpPos" structs snp_pos
    vector<SNPPOS> snp_pos;
    
    size_t pos;
    string rs;
    string chr;
    long int bp;
    size_t tt=0;
    vector<bool> indicator_func;
    vector<double> weight;
    double weight_i;
    
    for (size_t i=0; i < ns_total; ++i){
        if(indicator_snp[i] == 0) {continue;}
        
        if(tt == pos_loglr[tt].first ) pos = tt;
        else cout << "error assigning position to snp_pos vector"<< endl;
        
        rs = snpInfo[i].rs_number;
        chr = snpInfo[i].chr;
        bp = snpInfo[i].base_position;
        indicator_func = snpInfo[i].indicator_func;
        weight = snpInfo[i].weight;
        
        weight_i = snpInfo[i].weight_i;
        SNPPOS snp_temp={pos, rs, chr, bp, indicator_func, weight, weight_i};
        snp_pos.push_back(snp_temp);
        
        tt++;
    }
    snpInfo.clear();
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // sort snp by chr and bp
    //end of Jingjing's edit
    stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio
    
    //Jingjing's edit, create maps between rank and order
    for (size_t i=0; i<ns_test; ++i) {
        mapRank2pos[i]=pos_loglr[i].first;
        mapPos2Rank[pos_loglr[i].first] = i;
        
        mapOrder2pos[i] = snp_pos[i].pos;
        mapPos2Order[snp_pos[i].pos] = i;
    }
    
    for (size_t i=0; i<ns_test; ++i) {
        pos = mapRank2pos[i];
        mapRank2Order[i] = mapPos2Order[pos];
        
        pos = mapOrder2pos[i];
        mapOrder2Rank[i] = mapPos2Rank[pos];
    }
    //end of Jingjing's edit
    
    //calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
    gsl_rng_env_setup();
    const gsl_rng_type * gslType;
    gslType = gsl_rng_default;
    if (randseed<0)
    {
        time_t rawtime;
        time (&rawtime);
        tm * ptm = gmtime (&rawtime);
        
        randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
    }
    gsl_r = gsl_rng_alloc(gslType);
    gsl_rng_set(gsl_r, randseed);
    
    double *p_gamma = new double[ns_test];
    CalcPgamma (p_gamma); // calculate discrete distribution for gamma
    
    gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma); // set up proposal function for gamma
    
    //initial parameters
    cout << "Start initializing MCMC ... \n";
    InitialMCMC (X, z, rank_old, cHyp_old, pos_loglr, snp_pos); // Initialize rank and cHyp
    cHyp_initial=cHyp_old;
    CalcVPi(cHyp_old, pi_vec, sigma_vec, rank_old, snp_pos); // Calculate pi_vec, sigma_vec
    
    if (cHyp_old.n_gamma==0) {
        logPost_old=CalcPosterior (ztz, cHyp_old, pi_vec, sigma_vec, rank_old);
    }
    else {
        SetXgamma (Xgamma_old, X, rank_old);
        CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
        logPost_old=CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, pi_vec, sigma_vec, rank_old);
    }
    
    cout << "First 10 responses: \n";
    PrintVector(z, 10);
    cout << "First 10 Xtz: \n";
    PrintVector(Xtz_old, 10);
    cout << "First 10 X : \n";
    PrintMatrix(Xgamma_old, 10, 10);
    cout << "First 10 XtX : \n";
    PrintMatrix(XtX_old, 10, 10);
    
    //calculate centered z_hat, and pve
    if (a_mode==13) {
        if (cHyp_old.n_gamma==0) {
            CalcCC_PVEnZ (z_hat, cHyp_old);
        }
        else {
            CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
        }
    }
    
    //start MCMC
    int accept;
    size_t total_step=w_step+s_step;
    size_t w=0, w_col; // pos declared earlier (JY)
    //	size_t repeat=0;
    size_t repeat=1;
    int flag_gamma=0;
    double accept_percent;
    
    for (size_t t=0; t<total_step; ++t) {
        if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));}
        //		if (t>10) {break;}
        if (a_mode==13) {
            SampleZ (y, z_hat, z);
            mean_z=CenterVector (z);
            gsl_blas_ddot(z,z,&ztz);
            
            //First proposal
            if (cHyp_old.n_gamma==0) {
                logPost_old=CalcPosterior (ztz, cHyp_old, pi_vec, sigma_vec, rank_old);
            } else {
                gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
                gsl_vector_view Xtz_sub=gsl_vector_subvector(Xtz_old, 0, rank_old.size());
                gsl_blas_dgemv (CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0, &Xtz_sub.vector);
                logPost_old=CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, pi_vec, sigma_vec, rank_old);
            }
        }

        
        //MH steps
        //repeat = 1;
        //cout << "n_mh = " << n_mh << endl;
        for (size_t i=0; i<n_mh; ++i) {
            //if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
            //else {repeat=1;}
            
            logMHratio=0.0;
           // logMHratio+=ProposeHnRho(cHyp_old, cHyp_new, repeat);
           // logMHratio+=ProposePi(cHyp_old, cHyp_new, repeat);
          //  cout << "propose h, rho, pi success, proposing gamma..." << endl;
            
            logMHratio+=ProposeGamma (rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat, X, z, Xgamma_old, XtX_old, Xtz_old,  ztz, flag_gamma); //JY
            if(flag_gamma==1) nadd++;
            else if(flag_gamma==2) ndel++;
            else if(flag_gamma==3) nswitch++;
            else nother++;
          //  cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
            
            if (cHyp_old.n_gamma==0) {
                logPost_old=CalcPosterior (ztz, cHyp_old, pi_vec, sigma_vec, rank_old);
            }
            else {
                SetXgamma (Xgamma_old, X, rank_old);
                CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
                logPost_old=CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, pi_vec, sigma_vec, rank_old);
            }
            
            if (cHyp_new.n_gamma==0) {
                logPost_new=CalcPosterior (ztz, cHyp_new, pi_vec, sigma_vec, rank_new);
                
            } else {
                //this if makes sure that rank_old.size()==rank_remove.size() does not happen
                if (cHyp_new.n_gamma<=20 || cHyp_old.n_gamma<=20) {
                    time_start=clock();
                    SetXgamma (Xgamma_new, X, rank_new);
                    CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
                    time_set+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
                } else {
                    time_start=clock();
                    //cout << "start set Xgamma_new... " << endl;
                    SetXgamma (X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_new, XtX_new, Xtz_new);
                    time_set+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
                }
               // cout << "start calc posterior... " << endl;
                time_start=clock();
                logPost_new=CalcPosterior (Xgamma_new, XtX_new, Xtz_new, ztz, Xb_new, beta_new, cHyp_new, pi_vec, sigma_vec, rank_new);
                time_post+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            }
           //  cout << "Calcposterior success." << endl;
            logMHratio+=logPost_new-logPost_old;
           // cout << "logMHratio = " << logMHratio << endl;
            
            if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio) {accept=1; n_accept++;}
            else {accept=0;}
           // cout << "accept = " << accept << endl;
            
            if (accept==1) {
                
                if(flag_gamma==1) nadd_accept++;
                else if(flag_gamma==2) ndel_accept++;
                else if(flag_gamma==3) nswitch_accept++;
                else nother_accept++;
                
                //cout << "accept" << endl;
                logPost_old=logPost_new;
                cHyp_old=cHyp_new;
                gsl_vector_memcpy (Xb_old, Xb_new);
                
                rank_old.clear();
                if (rank_new.size()!=0) {
                    for (size_t ii=0; ii<rank_new.size(); ++ii) {
                        rank_old.push_back(rank_new[ii]);
                    }
                    
                    gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXold_sub=gsl_matrix_submatrix(XtX_old, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtzold_sub=gsl_vector_subvector(Xtz_old, 0, rank_new.size());
                    gsl_vector_view betaold_sub=gsl_vector_subvector(beta_old, 0, rank_new.size());
                    
                    gsl_matrix_view Xnew_sub=gsl_matrix_submatrix(Xgamma_new, 0, 0, ni_test, rank_new.size());
                    gsl_matrix_view XtXnew_sub=gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
                    gsl_vector_view Xtznew_sub=gsl_vector_subvector(Xtz_new, 0, rank_new.size());
                    gsl_vector_view betanew_sub=gsl_vector_subvector(beta_new, 0, rank_new.size());
                    
                    gsl_matrix_memcpy(&Xold_sub.matrix, &Xnew_sub.matrix);
                    gsl_matrix_memcpy(&XtXold_sub.matrix, &XtXnew_sub.matrix);
                    gsl_vector_memcpy(&Xtzold_sub.vector, &Xtznew_sub.vector);
                    gsl_vector_memcpy(&betaold_sub.vector, &betanew_sub.vector);
                }
            } else {
                cHyp_new=cHyp_old;
            }
            //cout << "copy data from new propose -> old " << endl;
        }
        
        //calculate z_hat, and pve
        if (a_mode==13) {
            if (cHyp_old.n_gamma==0) {
                CalcCC_PVEnZ (z_hat, cHyp_old);
            }
            else {
                CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
            }
            
            //sample mu and update z hat
            gsl_vector_sub (z, z_hat);
            mean_z+=CenterVector(z);
            mean_z+=gsl_ran_gaussian(gsl_r, sqrt(1.0/(double) ni_test) );
            
            gsl_vector_add_constant (z_hat, mean_z);
        }
        
        accept_percent = (double)n_accept/(double)(t*n_mh);
        if (accept_percent<0.000001) {
            cerr << "acceptance percentage = " << accept_percent << " < 0.000001; ABORT MCMC...";
            exit(1);
        }
        
        //Save data
        if (t<w_step) {continue;}
        else {
            if (t%r_pace==0) {
                w_col=w%w_pace;
                if (w_col==0) {
                    if (w==0) {
                        WriteResult (0, Result_hyp, Result_gamma, Result_theta, Result_sigma, w_col);
                       // WriteResult (0, Result_hyp, Result_gamma, w_col, snp_pos, pos_loglr);
                    }
                    else {
                        WriteResult (1, Result_hyp, Result_gamma, Result_theta, Result_sigma,w_col);
                       // WriteResult (1, Result_hyp, Result_gamma, w_col, snp_pos, pos_loglr);
                        gsl_matrix_set_zero (Result_hyp);
                        gsl_matrix_set_zero (Result_gamma);
                        gsl_matrix_set_zero (Result_theta);
                        gsl_matrix_set_zero (Result_sigma);
                    }
                }
                
                gsl_matrix_set (Result_hyp, w_col, 0, cHyp_old.pve);
                gsl_matrix_set (Result_hyp, w_col, 1, cHyp_old.pge);
                gsl_matrix_set (Result_hyp, w_col, 2, cHyp_old.n_gamma);
                gsl_matrix_set (Result_hyp, w_col, 3, logPost_old);
                
                for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
                    // beta_g saved by position
                    
                    pos=mapRank2pos[rank_old[i]]+1;
                    
                    gsl_matrix_set (Result_gamma, w_col, i, pos);
                    
                    beta_g[pos-1].first+=gsl_vector_get(beta_old, i);
                    beta_g[pos-1].second+=1.0;	
                }
                
                if (a_mode==13) {
                    pheno_mean+=mean_z;
                }
                
                w++;
                
            }
            
        }
    }
    cout<<endl;
    
    cout<<"time on selecting Xgamma: "<<time_set<<endl;
    cout<<"time on calculating posterior: "<<time_post<<endl;
    
    w_col=w%w_pace;
   // WriteResult (1, Result_hyp, Result_gamma, w_col, snp_pos, pos_loglr);
    WriteResult (1, Result_hyp, Result_gamma, Result_theta, Result_sigma, w_col);
    
    gsl_vector *alpha=gsl_vector_alloc (ns_test);
    gsl_vector_set_zero (alpha);
    //WriteParam (beta_g, alpha, w);
    WriteParam (beta_g, alpha, w, snp_pos, pos_loglr);
    gsl_vector_free(alpha);
    
    gsl_matrix_free(Result_hyp);
    gsl_matrix_free(Result_gamma);
    gsl_matrix_free(Result_theta);
    gsl_matrix_free(Result_sigma);
    
    gsl_vector_free(pi_vec);
    gsl_vector_free(sigma_vec);
    
    gsl_vector_free(z_hat);
    gsl_vector_free(z);
    gsl_vector_free(Xb_new);	
    gsl_vector_free(Xb_old);
    
    gsl_matrix_free(Xgamma_old);
    gsl_matrix_free(XtX_old);
    gsl_vector_free(Xtz_old);
    gsl_vector_free(beta_old);
    
    gsl_matrix_free(Xgamma_new);
    gsl_matrix_free(XtX_new);
    gsl_vector_free(Xtz_new);
    gsl_vector_free(beta_new);
    
    //gsl_ran_discrete_free(gsl_t);
    
    delete [] p_gamma;
    beta_g.clear();
    
    return;
}
// end from previous version


//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	
	gsl_blas_ddot (Utu, Utu, &d);
	cHyp.pve=d/(double)ni_test;	
		
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu, 0.0, z_hat);
		
	cHyp.pve/=cHyp.pve+1.0;
	cHyp.pge=0.0;	
	
	return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *UtXb, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	gsl_vector *UtXbU=gsl_vector_alloc (Utu->size);
	
	gsl_blas_ddot (UtXb, UtXb, &d);
	cHyp.pge=d/(double)ni_test;
	
	gsl_blas_ddot (Utu, Utu, &d);
	cHyp.pve=cHyp.pge+d/(double)ni_test;
	
	gsl_vector_memcpy (UtXbU, Utu);
	gsl_vector_add (UtXbU, UtXb);
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, UtXbU, 0.0, z_hat);	
	
	if (cHyp.pve==0) {cHyp.pge=0.0;}
	else {cHyp.pge/=cHyp.pve;}
	
	cHyp.pve/=cHyp.pve+1.0;
	
	gsl_vector_free(UtXbU);
	return;
}


void BSLMM::SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z)
{	
	double d1, d2, z_rand=0.0;
	for (size_t i=0; i<z->size; ++i) {
		d1=gsl_vector_get (y, i);
		d2=gsl_vector_get (z_hat, i);
		//y is centerred for case control studies
		if (d1<=0.0) {
			//control, right truncated
			do {				
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand>0.0);
		}
		else {
			do {
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand<0.0);
		}
		
		gsl_vector_set (z, i, z_rand);
	}

	return;
}


double BSLMM::ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	
	double h=cHyp_old.h, rho=cHyp_old.rho;
	
	double d_h=(h_max-h_min)*h_scale, d_rho=(rho_max-rho_min)*rho_scale;
	
	for (size_t i=0; i<repeat; ++i) {
		h=h+(gsl_rng_uniform(gsl_r)-0.5)*d_h;
		if (h<h_min) {h=2*h_min-h;}
		if (h>h_max) {h=2*h_max-h;}
		
		rho=rho+(gsl_rng_uniform(gsl_r)-0.5)*d_rho;
		if (rho<rho_min) {rho=2*rho_min-rho;}
		if (rho>rho_max) {rho=2*rho_max-rho;}
	}
	/*
	//Grid Sampling
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		h=h+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.1;
		if (h<h_min) {h=h_max;}
		if (h>h_max) {h=h_min;}
	}
	
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		rho=rho+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.1;
		if (rho<rho_min) {rho=rho_max;}
		if (rho>rho_max) {rho=rho_min;}
	}
	*/
	cHyp_new.h=h;
	cHyp_new.rho=rho;
	return 0.0;
}


double BSLMM::ProposePi (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	double logp_old=cHyp_old.logp, logp_new=cHyp_old.logp;
	double log_ratio = 0.0;
    
	double d_logp=min(0.1, (logp_max-logp_min)*logp_scale);
	
	for (size_t i=0; i<repeat; ++i) {
        
		logp_new=logp_old+(gsl_rng_uniform(gsl_r)-0.5)*d_logp;
		if (logp_new<logp_min) {logp_new=2*logp_min-logp_new;}
		if (logp_new>logp_max) {logp_new=2*logp_max-logp_new;}		
		
		log_ratio+=logp_new-logp_old;
		logp_old=logp_new;
        
	}
	/*
	//Grid Sampling
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		logp_new=logp_old+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.5*log(10.0);
		if (logp_new<logp_min) {logp_new=logp_max;}
		if (logp_new>logp_max) {logp_new=logp_min;}	
		
		log_ratio+=logp_new-logp_old;
		logp_old=logp_new;
	}
	*/
	cHyp_new.logp=logp_new;
	
	return log_ratio;
}

double BSLMM::ProposeTheta (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	vector<double> logp_old = cHyp_old.log_theta, logp_new=cHyp_old.log_theta;
    
	double log_ratio = 0.0;
    
	double d_logp = min(0.1, (logp_max-logp_min)*logp_scale);
	
	for (size_t i=0; i<repeat; ++i) {
        
        for (size_t j=0; j<n_type; j++) {
   
            logp_new[j] = logp_old[j] + (gsl_rng_uniform(gsl_r)-0.5)*d_logp;
            if (logp_new[j]<logp_min) {logp_new[j]=2*logp_min-logp_new;}
            if (logp_new[j]>logp_max) {logp_new[j]=2*logp_max-logp_new;}
            log_ratio+=logp_new[j]-logp_old[j];
            logp_old[j]=logp_new[j];
        }
	}
	cHyp_new.log_theta=logp_new;
	
	return log_ratio;
}

double BSLMM::ProposeSubvar (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	vector<double> logp_old = cHyp_old.subvar, logp_new=cHyp_old.subvar;
    
	double log_ratio = 0.0;
    	
	for (size_t i=0; i<repeat; ++i) {
        
        for (size_t j=0; j<n_type; j++) {
            
            logp_new[j] = logp_old[j] + gsl_ran_gaussian_tail(gsl_r, -logp_old[j], sigma_prop);
            log_ratio += log(gsl_cdf_gaussian_Q(-logp_old[j], sigma_prop)) - log(gsl_cdf_gaussian_Q(-logp_new[j], sigma_prop));
            logp_old[j]=logp_new[j];
        }
	}
	cHyp_new.subvar=logp_new;
	
	return log_ratio;
}


//JY edit start

void BSLMM::CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX, const gsl_vector *Xtz, gsl_vector *z_res, const size_t &s_size, const double &ztz){
    
    gsl_matrix_const_view X_gsub=gsl_matrix_const_submatrix(Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_gsub = gsl_matrix_const_submatrix(XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xtz_gsub = gsl_vector_const_subvector(Xtz, 0, s_size);
    
    gsl_vector *beta_gamma_hat = gsl_vector_alloc(s_size);
    
    double lambda = 0.0;
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(XtX, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.00000001;
   // cout << "labmda = " << lambda << endl;
    
    EigenSolve(&XtX_gsub.matrix, &Xtz_gsub.vector, beta_gamma_hat, lambda);
    //EigenSolve(&XtX_gsub.matrix, &Xtz_gsub.vector, beta_gamma_hat);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &X_gsub.matrix, beta_gamma_hat, 0.0, z_res);
    gsl_vector_scale(z_res, -1.0);
    gsl_vector_add(z_res, z);
    
    double SSR;
    gsl_blas_ddot(z_res, z_res, &SSR);
    double R2 = 1.0 - (SSR / ztz);
    // R2 = R2 - (1 - R2) * s_size / (ni_test - s_size - 1);
   // cout << "R2 = "<< R2 << endl;
    if(R2 <= 0.0) {
        cout << "R2 = " << setprecision(6) << R2 << " <= 0, ";
        cout << "beta_hat estimate from calculating residuals: \n ";
        PrintVector(beta_gamma_hat);
        //cout << "SSR = " << SSR << "; ztz = " << ztz << "\n";
        //cout << "Set z_res equal to z... " << endl;
        gsl_vector_memcpy(z_res, z);
        //NormRes(z_res);
        }
        
    double res_mean = CenterVector(z_res);
    if(res_mean > 0.0001) cout << "res_mean = " << res_mean  << " > 0.0001" << endl;
    
    gsl_vector_free(beta_gamma_hat);
    return;
}

bool comp_res (pair<size_t, double> a, pair<size_t, double> b)
{
	return (a.second < b.second);
}

void BSLMM::NormRes(gsl_vector * z_res){
    
    size_t vec_length = (z_res->size);
    size_t y_ind;
    vector<pair<size_t, double> > y_res;
    //vector<double> y_norm;
    
    for (size_t i=0; i<vec_length; ++i) {
        y_res.push_back(make_pair(i, gsl_vector_get(z_res, i)));
    }
    std::random_shuffle(y_res.begin(), y_res.end());
    stable_sort (y_res.begin(), y_res.end(), comp_res);
    
    double unit_p = 1.0/(vec_length+1);
    double p_i = unit_p, qnorm;
    for (size_t i=0; i<vec_length; ++i) {
        qnorm = gsl_cdf_ugaussian_Pinv(p_i);
        p_i += unit_p;
        y_ind = y_res[i].first;
        gsl_vector_set(z_res, y_ind, qnorm);
    }
   // cout << "normalize success" << endl;
    return;
}

double BSLMM::CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec){
    double LR;
    double xtx_vec, xtz_res, ztz_res;
    
    gsl_blas_ddot(z_res, z_res, &ztz_res);
    gsl_blas_ddot(x_vec, x_vec, &xtx_vec);
    gsl_blas_ddot(x_vec, z_res, &xtz_res);
    //cout << "ztz_res = " << ztz_res << "; xtx_vec = " << xtx_vec << "; xtz_res = " << xtz_res << endl;
    //double ixtx = 1.0 / xtx_vec;
    //double bhat = ixtx * xtz_res;
    //double V = (ixtx * ztz_res - bhat * bhat) / (ni_test);
    //double Z2 = bhat * bhat / V;
    //double VpW = V + Wvar;
    LR = 0.5*(ni_test)*(log(ztz_res)-log(ztz_res-xtz_res*xtz_res/xtx_vec));
    //sqrt(VpW / V) * exp(-0.5 * Z2 * Wvar / VpW);
    //cout << "log LR = " << BF << ", ";
    return (LR);
}

gsl_ran_discrete_t * BSLMM::MakeProposal(const size_t &o, double *p_BF, uchar **X, const gsl_vector *z_res, const map<size_t, int> &mapRank2in)
{
    gsl_vector *xvec = gsl_vector_alloc(z_res->size);
    
    long int orderj;
    size_t posj, rank_j;
    vector<int> j_ind;
    double pmax, psum=0.0, countj = 0.0;
    
    for (size_t j=0; j < ns_neib; ++j){
        orderj = (o - win) + j;
        rank_j = mapOrder2Rank[orderj];
        if((orderj >= 0) && (orderj < (long int)ns_test) && (j != win) && (mapRank2in.count(rank_j) == 0)){
            posj = mapOrder2pos[orderj];
            getGTgslVec(X, xvec, posj, ni_test, ns_test, CompBuffSizeVec, UnCompBufferSize);
            p_BF[j]=CalcLR(z_res, xvec); //calc loglr
            j_ind.push_back(1);
            countj += 1.0;
        }
        else{
            p_BF[j] = -std::numeric_limits<double>::max();
            j_ind.push_back(0);
        }
    }
    pmax = *std::max_element(p_BF, p_BF+(ns_neib));
    
    for (size_t j=0; j < ns_neib; ++j){
        orderj = (o - win) + j;
        if(j_ind[j]==1){
            p_BF[j]=exp(p_BF[j]- pmax);
            psum += p_BF[j];
        }
        else{p_BF[j] = 0.0;}
    }
    //
    psum = 1.0/psum;
    for(size_t j=0; j < ns_neib; ++j){
        p_BF[j] *= psum;
    }
    
    gsl_vector_free(xvec);
    
    return (gsl_ran_discrete_preproc(ns_neib, p_BF));
}


// JY edit end

void BSLMM::AddMarker(double &logp, map<size_t, int> &mapRank2in, class HYPBSLMM &cHyp_new, vector<size_t> &rank_new, const double *p_gamma){
    
    size_t r_add, r;

    do {
        r_add=gsl_ran_discrete (gsl_r, gsl_t);
    } while (mapRank2in.count(r_add)!=0);
    
    double prob_total=1.0;
    for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
        r=rank_new[i];
        prob_total-=p_gamma[r];
    }
    
    mapRank2in[r_add]=1;
    rank_new.push_back(r_add);
    cHyp_new.n_gamma++;
    logp+=-log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);
    
    return;
}

void BSLMM::DelMarker(double &logp, map<size_t, int> &mapRank2in, class HYPBSLMM &cHyp_new, vector<size_t> &rank_new, const double *p_gamma)
{
    size_t r_remove, col_id, r;

    col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
    r_remove=rank_new[col_id];
    
    double prob_total=1.0;
    for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
        r=rank_new[i];
        prob_total-=p_gamma[r];
    }
    prob_total+=p_gamma[r_remove];
    
    mapRank2in.erase(r_remove);
    rank_new.erase(rank_new.begin()+col_id);
    logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
    cHyp_new.n_gamma--;
    return;
}

void BSLMM::SwitchMarker(double &logp, map<size_t, int> &mapRank2in, LModel &model_old, LModel &model_new, uchar **X, gsl_vector *z, const double &ztz){
    
    size_t r_add, r_remove, col_id, r;
    long int o_add, o_remove;
    long int o_rj, o_aj;
    size_t j_add, j_remove, o;
    
    LModel model_temp; // used for making proposal distribution
    model_temp.InitialVar(ni_test, s_max);
    
    gsl_vector *z_res = gsl_vector_alloc(ni_test);
    gsl_ran_discrete_t *gsl_s, *gsl_a; //JY added dynamic gsl_s
    
    double *p_BFr = new double[ns_neib];
    double *p_BFa = new double[ns_neib];
    
    //remove a marker
    col_id=gsl_rng_uniform_int(gsl_r, model_new.cHyp.n_gamma);
    r_remove=model_new.rank[col_id];//careful with the proposal
    if(mapRank2in.count(r_remove) == 0) {cout << "wrong proposal of r_remove;" << endl; exit(1);}
    o_remove = mapRank2Order[r_remove];
    model_new.rank.erase(model_new.rank.begin()+col_id);
    size_t s_size = model_new.rank.size();
    mapRank2in.erase(r_remove);
    
    
    model_temp.Copy(model_new);
    //model_temp.cHyp.n_gamma = model_new.cHyp.n_gamma;
    if (model_temp.cHyp.n_gamma<=20 || model_old.cHyp.n_gamma<=20) {
        model_temp.AssignVar(X, z, mapRank2pos, ns_test);
    } else {
        SetXgamma (model_old, model_temp, X, z);
    }
    CalcRes(model_temp.Xgamma, z, model_temp.XtX, model_temp.Xty, z_res, s_size, ztz);
    gsl_s = MakeProposal(o_remove, p_BFr, X, z_res, mapRank2in);
    //construct gsl_s, JY
    
    j_add = gsl_ran_discrete(gsl_r, gsl_s);
    o_add = (o_remove - win) + j_add;
    r_add = mapOrder2Rank[o_add];
    //cout << "o_add = " << o_add <<  "; r_add = "<<r_add << endl;
    if((o_add < 0) || (o_add >= (long int)ns_test) || (o_add == (long int)o_remove))
        cout << "ERROR proposing switch snp"; //new snp != removed snp
    
    gsl_a = MakeProposal(o_add, p_BFa, X, z_res, mapRank2in);
    
    double prob_total_remove=1.0;
    double prob_total_add=1.0;
    for (size_t i=0; i<model_new.rank.size(); ++i) {
        r = model_new.rank[i];
        o = mapRank2Order[r];
        o_rj = ((long int)o - o_remove) + win;
        o_aj = ((long int)o - o_add) + win;
        if(o_aj >= 0 && o_aj < (long int)ns_neib) prob_total_add -= p_BFa[o_aj];
        if(o_rj >= 0 && o_rj < (long int)ns_neib) prob_total_remove -= p_BFr[o_rj];
    }
    j_remove = o_remove - o_add + win;
    logp += log( p_BFa[j_remove] / prob_total_add ); //prob(delete o_add & add o_remove)
    logp -= log( p_BFr[j_add] / prob_total_remove ); //prob(delete o_remove & add o_add)
    
    mapRank2in[r_add]=1;
    model_new.rank.push_back(r_add);
    
    model_temp.FreeMem(); // free model_temp
    gsl_vector_free(z_res);
    gsl_ran_discrete_free(gsl_s);
    gsl_ran_discrete_free(gsl_a);
    delete[] p_BFr;
    delete[] p_BFa;
    
    return;

}

void BSLMM::AssignRank(map<size_t, int> &mapRank2in, vector<size_t> &rank_new, const vector<size_t> &rank_old, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new)
{
    size_t r;
    rank_new.clear();
    
	if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"size wrong"<<endl;}
	if (cHyp_old.n_gamma!=0) {
		for (size_t i=0; i<rank_old.size(); ++i) {
			r=rank_old[i];
			rank_new.push_back(r);
			mapRank2in[r]=1;
		}
	}
    
    cHyp_new.n_gamma=cHyp_old.n_gamma;
    return;
}


double BSLMM::ProposeGamma (LModel &model_old, LModel &model_new, const double *p_gamma, const size_t &repeat, uchar **X, gsl_vector *z, const double &ztz, int &flag_gamma)
{
	map<size_t, int> mapRank2in;
	double unif, logp = 0.0;
    
    //set rank_new = rank_old; cHyp_new.n_gamma=cHyp_old.n_gamma;
	AssignRank(mapRank2in, model_new.rank, model_old.rank, model_old.cHyp, model_new.cHyp);
	    
	for (size_t i=0; i<repeat; ++i) {
        
		unif=gsl_rng_uniform(gsl_r);
        
		if (unif < 0.33 && model_new.cHyp.n_gamma<s_max)
        {
            flag_gamma=1;
            AddMarker(logp, mapRank2in, model_new.cHyp, model_new.rank, p_gamma);
        }
		else if (unif>=0.33 && unif < 0.67 && model_new.cHyp.n_gamma>s_min)
        {
            flag_gamma=2;
            DelMarker(logp, mapRank2in, model_new.cHyp, model_new.rank, p_gamma);
        }
		else if (unif>=0.67 && model_new.cHyp.n_gamma>0 && model_new.cHyp.n_gamma<ns_test)
        {
            flag_gamma=3;
            SwitchMarker(logp, mapRank2in, model_old, model_new, X, z, ztz);
        }
		else {flag_gamma=4;}
	}
	
	stable_sort (model_new.rank.begin(), model_new.rank.end(), comp_vec);
	mapRank2in.clear();
	return logp;
    
}




void BSLMM::RidgeR(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *eval, const double lambda)
{
	gsl_vector *beta=gsl_vector_alloc (UtX->size2);
	gsl_vector *H_eval=gsl_vector_alloc (Uty->size);
	gsl_vector *bv=gsl_vector_alloc (Uty->size);

	gsl_vector_memcpy (H_eval, eval);
	gsl_vector_scale (H_eval, lambda);
	gsl_vector_add_constant (H_eval, 1.0);
	
	gsl_vector_memcpy (bv, Uty);
	gsl_vector_div (bv, H_eval);	

	gsl_blas_dgemv (CblasTrans, lambda/(double)UtX->size2, UtX, bv, 0.0, beta);
	gsl_vector_add_constant (H_eval, -1.0);
	gsl_vector_mul (H_eval, bv);
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, H_eval, 0.0, bv);

	WriteParam (beta);
	WriteBV(bv);
	
	gsl_vector_free (H_eval);
	gsl_vector_free (beta);
	gsl_vector_free (bv);
	
	return;
}



//below fits MCMC for rho=1
void BSLMM::CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty)
{
  time_t time_start=clock();	
  gsl_matrix_const_view X_sub=gsl_matrix_const_submatrix(X, 0, 0, X->size1, s_size);
  gsl_matrix_view XtX_sub=gsl_matrix_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_view Xty_sub=gsl_vector_subvector(Xty, 0, s_size);

  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &X_sub.matrix, &X_sub.matrix, 0.0, &XtX_sub.matrix);
  gsl_blas_dgemv(CblasTrans, 1.0, &X_sub.matrix, y, 0.0, &Xty_sub.vector);

  time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

  return;
}


void BSLMM::SetXgamma (LModel &model_old, LModel &model_new, uchar **X, const gsl_vector *y)
{
  double d;
   // cout << "X_add set start" << endl;

  //rank_old and rank_new are sorted already inside PorposeGamma
  //calculate vectors rank_remove and rank_add
  //  size_t v_size=max(rank_old.size(), rank_new.size());
  //make sure that v_size is larger than repeat
  size_t v_size=20;
  vector<size_t> rank_remove(v_size), rank_add(v_size), rank_union(s_max+v_size);
  vector<size_t>::iterator it;

  it=set_difference (model_old.rank.begin(), model_old.rank.end(), model_new.rank.begin(), model_new.rank.end(), rank_remove.begin());
  rank_remove.resize(it-rank_remove.begin());

  it=set_difference (model_new.rank.begin(), model_new.rank.end(), model_old.rank.begin(), model_old.rank.end(), rank_add.begin());
  rank_add.resize(it-rank_add.begin());

  it=set_union (model_new.rank.begin(), model_new.rank.end(), model_old.rank.begin(), model_old.rank.end(), rank_union.begin());
  rank_union.resize(it-rank_union.begin());

  //map rank_remove and rank_add
  map<size_t, int> mapRank2in_remove, mapRank2in_add;
  for (size_t i=0; i<rank_remove.size(); i++) {
    mapRank2in_remove[rank_remove[i]]=1;
  }
  for (size_t i=0; i<rank_add.size(); i++) {
    mapRank2in_add[rank_add[i]]=1;
  }

    
    model_old.getSubVar(model_old.rank.size());
    model_new.getSubVar(model_new.rank.size());
  //obtain the subset of matrix/vector

   if (rank_remove.size()==0 && rank_add.size()==0) {
       model_new.Copy(model_old);
     // cout << "rank_old = rank_new; " << "Xgamma_new set success" << endl;
  } else {
    size_t i_old, j_old, i_new, j_new, i_add, j_add, i_flag, j_flag;
    if (rank_add.size()==0) {
      i_old=0; i_new=0;
      for (size_t i=0; i<rank_union.size(); i++) {
	if (mapRank2in_remove.count(model_old.rank[i_old])!=0) {i_old++; continue;}

	gsl_vector_view Xnew_col=gsl_matrix_column(model_new.Xgamma, i_new);
	gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(model_old.Xgamma, i_old);
	gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);

	d=gsl_vector_get (model_old.Xty, i_old);
	gsl_vector_set (model_new.Xty, i_new, d);

	j_old=i_old; j_new=i_new;
	for (size_t j=i; j<rank_union.size(); j++) {
      if (mapRank2in_remove.count(model_old.rank[j_old])!=0) {j_old++; continue;}

	  d=gsl_matrix_get(model_old.XtX, i_old, j_old);
	  gsl_matrix_set (model_new.XtX, i_new, j_new, d);
	  if (i_new!=j_new) {gsl_matrix_set (model_new.XtX, j_new, i_new, d);}
            j_old++; j_new++;
        }
          i_old++; i_new++;
      }
       // cout << "X_add = NULL; " << "Xgamma_new set success" << endl;
    } else {
        //rank_add has length > 0
      gsl_matrix *X_add=gsl_matrix_alloc(ni_test, rank_add.size() );
      gsl_matrix *XtX_aa=gsl_matrix_alloc(X_add->size2, X_add->size2);
      gsl_matrix *XtX_ao=gsl_matrix_alloc(X_add->size2, model_old.rank.size());
      gsl_vector *Xty_add=gsl_vector_alloc(X_add->size2);

      //get X_add
      SetXgamma (X_add, X, rank_add);
      //get t(X_add)X_add and t(X_add)X_temp	
      clock_t time_start=clock();
      gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, X_add, 0.0, XtX_aa);
      gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X_add, &model_old.Xgamma_sub.matrix, 0.0, XtX_ao);
      gsl_blas_dgemv(CblasTrans, 1.0, X_add, y, 0.0, Xty_add);
      time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

      //save to X_new, XtX_new and Xty_new
      i_old=0; i_new=0; i_add=0;
      for (size_t i=0; i<rank_union.size(); i++) {
	if (mapRank2in_remove.count(model_old.rank[i_old])!=0) {i_old++; continue;}
	if (mapRank2in_add.count(model_new.rank[i_new])!=0){i_flag=1;} //within x_add
    else {i_flag=0;}//within x_common

	gsl_vector_view Xnew_col=gsl_matrix_column(model_new.Xgamma, i_new);
	if (i_flag==1) {
	  gsl_vector_view Xcopy_col=gsl_matrix_column(X_add, i_add);
	  gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
	} else {
	  gsl_vector_const_view Xcopy_col=gsl_matrix_const_column(model_old.Xgamma, i_old);
	  gsl_vector_memcpy (&Xnew_col.vector, &Xcopy_col.vector);
	}

	if (i_flag==1) {
          d=gsl_vector_get (Xty_add, i_add);
        } else {
          d=gsl_vector_get (model_old.Xty, i_old);
        }          
	gsl_vector_set (model_old.Xty, i_new, d);

	j_old=i_old; j_new=i_new; j_add=i_add;
	for (size_t j=i; j<rank_union.size(); j++) {
	  if (mapRank2in_remove.count(model_old.rank[j_old])!=0) {j_old++; continue;} // within rank_remove
	  if (mapRank2in_add.count(model_new.rank[j_new])!=0) {j_flag=1;}//within rank_add
      else {j_flag=0;} // within rank_old

	  if (i_flag==1 && j_flag==1) {
            d=gsl_matrix_get(XtX_aa, i_add, j_add);          
	  } else if (i_flag==1) {
	    d=gsl_matrix_get(XtX_ao, i_add, j_old);
	  } else if (j_flag==1) {
	    d=gsl_matrix_get(XtX_ao, j_add, i_old);
	  } else {
	    d=gsl_matrix_get(model_old.XtX, i_old, j_old);
	  }

	  gsl_matrix_set (model_new.XtX, i_new, j_new, d);
	  if (i_new!=j_new) {gsl_matrix_set (model_new.XtX, j_new, i_new, d);}

	  j_new++; if (j_flag==1) {j_add++;} else {j_old++;}
        }
          //cout << "XtX_new success" << endl;
	i_new++; if (i_flag==1) {i_add++;} else {i_old++;}
      }
       // cout << "X_gamma set success" << endl;
        
      gsl_matrix_free(X_add);
      gsl_matrix_free(XtX_aa);
      gsl_matrix_free(XtX_ao);
      gsl_vector_free(Xty_add);
    }

  }

  rank_remove.clear();
  rank_add.clear();
  rank_union.clear();
  mapRank2in_remove.clear();
  mapRank2in_add.clear();
	
  return;
}


void BSLMM::LModel::InitialVar(size_t &ni_test, size_t &s_max)
{
    Xgamma = gsl_matrix_alloc (ni_test, s_max);
    XtX = gsl_matrix_alloc (s_max, s_max);
    Xty = gsl_vector_alloc (s_max);
    beta = gsl_vector_alloc (s_max);
    Xbeta=gsl_vector_alloc (ni_test);
    return;
}

void BSLMM::LModel::getSubVar(size_t subsize)
{
    Xgamma_sub=gsl_matrix_submatrix (Xgamma, 0, 0, Xgamma->size1, subsize);
    XtX_sub=gsl_matrix_submatrix (XtX, 0, 0, subsize, subsize);
    Xty_sub=gsl_vector_subvector (Xty, 0, subsize);
    beta_sub=gsl_vector_subvector (beta, 0, subsize);
}

void BSLMM::LModel::SetXgamma (uchar **X, map<size_t, size_t> &mapRank2pos, size_t &ns_test)
{
    size_t pos;
    for (size_t i=0; i<rank.size(); ++i) {
        pos=mapRank2pos[rank[i]];
        gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
        getGTgslVec(X, &Xgamma_col.vector, pos, Xgamma->size1, ns_test);
    }
    return;
}

void BSLMM::LModel::AssignVar(uchar **X, gsl_vector *y, map<size_t, size_t> &mapRank2pos, size_t &ns_test)
{
    SetXgamma(X, mapRank2pos, ns_test);
    getSubVar(rank.size());
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &Xgamma_sub.matrix, &Xgamma_sub.matrix, 0.0, &XtX_sub.matrix);
    gsl_blas_dgemv(CblasTrans, 1.0, &Xgamma_sub.matrix, y, 0.0, &Xty_sub.vector);
    
}

void BSLMM::LModel::Copy(LModel &model)
{
    cHyp = model.cHyp;
    gsl_vector_memcpy(Xbeta, model.Xbeta);
    
    rank.clear();
    if (model.rank.size()!=0) {
        for (size_t i=0; i<model.rank.size(); ++i) {
            rank.push_back(model.rank[i]);
        }
        
        getSubVar(model.rank.size());
        model.getSubVar(model.rank.size());
        
        gsl_matrix_memcpy(&Xgamma_sub.matrix, &model.Xgamma_sub.matrix);
        gsl_matrix_memcpy(&XtX_sub.matrix, &model.XtX_sub.matrix);
        gsl_vector_memcpy(&Xty_sub.vector, &model.Xty_sub.vector);
        gsl_vector_memcpy(&beta_sub.vector, &model.beta_sub.vector);
    }
}

void BSLMM::LModel::FreeMem()
{
    gsl_matrix_free(Xgamma);
    gsl_matrix_free(XtX);
    gsl_vector_free(Xty);
    gsl_vector_free(beta);
    gsl_vector_free(Xbeta);
}

double BSLMM::LModel::CalcPosterior (const double yty, size_t &ni_test, size_t &ns_test, gsl_rng *gsl_r, int &a_mode, const double &trace_G)
{
    double logpost=0.0;
    
 if (rank.size()==0) {
        if (a_mode==11) {
            cHyp.pve=0.0;
            cHyp.pge=1.0;
        }
        //calculate likelihood
        if (a_mode==11) {logpost-=0.5*(double)ni_test*log(yty);}
        else {logpost-=0.5*yty;}
        
        logpost+=((double)(cHyp.n_gamma-1))*cHyp.logp+((double)(ns_test-cHyp.n_gamma))*log(1-exp(cHyp.logp));
    }
    
 else {
     
    // cout << "calculate postrior with ranks:\n";
     //PrintVector(rank);
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = rank.size();
     
    // cout << "h = " << cHyp.h << "; pi = " << exp(cHyp.logp) << "\n";
    // double sigma_a2=cHyp.h/((1.0-cHyp.h) * exp(cHyp.logp) * ((double)ns_test) * trace_G);
     double sigma_a2=cHyp.h/((1.0-cHyp.h) * exp(cHyp.logp) * ((double)ns_test));
     getSubVar(s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);
    gsl_vector *Xty_temp=gsl_vector_alloc (s_size);
    gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
    
    //calculate Omega
     gsl_matrix_memcpy (Omega, &XtX_sub.matrix);
    gsl_matrix_scale (Omega, sigma_a2);
    gsl_matrix_set_identity (M_temp);
    gsl_matrix_add (Omega, M_temp);
    
   //  cout << "Print out Omega sub-matrix:\n";
   //  PrintMatrix(Omega, 3, 10);
     
    //calculate beta_hat
   //  cout << "Print estimated Xty_temp : \n";
   //  PrintVector(Xty_temp) ;
     
    // EigenSolve(Omega, Xty_temp, beta_hat);
   //  cout << "Print estimated beta_hat from EigenSolve : \n";
   //  PrintVector(beta_hat);
     
    logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);	//solve Omega * beta_hat = Xty forbeta_hat
   //  cout << "Print estimated beta_hat from CholeskySolve: \n";
   //  PrintVector(beta_hat) ;
     
   //  cout << "logdet_0 = " << logdet_O << "\n";
    
    gsl_vector_scale (beta_hat, sigma_a2);
    // cout << "Print estimated beta_hat : \n";
    // PrintVector(beta_hat) ;
     
    gsl_blas_ddot (Xty_temp, beta_hat, &d);
    P_yy = yty - d;
    // cout << "; P_yy = " << P_yy << endl;
     if(P_yy <= 0) {cout << "Error: P_yy <= 0\n";
         cout << "cHyp.h = " << cHyp.h << "; exp(cHyp.logp) = " << exp(cHyp.logp);
         cout << "; sigma_a2=" << sigma_a2 << endl;
         
         EigenSolve(Omega, Xty_temp, beta_hat, 0.0000001);
         gsl_vector_scale (beta_hat, sigma_a2);
         gsl_blas_ddot (Xty_temp, beta_hat, &d);
         P_yy=yty - d;
     }
    
    //sample tau
    double tau=1.0;
    if (a_mode==11) {tau = gsl_ran_gamma (gsl_r, ((double)ni_test)/2.0,  2.0/P_yy); }
    // cout << "; tau = " << tau;
    
    //sample beta
    for (size_t i=0; i<s_size; i++)
    {
        d=gsl_ran_gaussian(gsl_r, 1);
        gsl_vector_set(beta, i, d);
    }
     
    //cout << "\n Sample beta =  \n";
    // PrintVector(&beta_sub.vector);
    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, &beta_sub.vector);
    
    gsl_vector_scale(&beta_sub.vector, sqrt(sigma_a2/tau));
    gsl_vector_add(&beta_sub.vector, beta_hat);
    gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xbeta);
    
    //for quantitative traits, calculate pve and pge
    if (a_mode==11) {
        gsl_blas_ddot (Xbeta, Xbeta, &d);
        cHyp.pve=d/(double)ni_test;
        cHyp.pve/=cHyp.pve+1.0/tau;
        cHyp.pge=1.0;
    }
     

    logpost= -0.5*logdet_O;
     //cout << "after logdet_0 logpost= " << logpost << "\n";

    if (a_mode==11) {logpost-=0.5*((double)ni_test)*log(P_yy);}
    else {logpost-=0.5*P_yy;}
    
    // cout << "after P_yy logpost= " << logpost << "\n";

    logpost+=((double)(cHyp.n_gamma-1))*cHyp.logp+((double)(ns_test-cHyp.n_gamma))*log(1.0-exp(cHyp.logp));
    
    gsl_matrix_free (Omega);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (Xty_temp);
    }
    // cout << "log posterior = " << logpost << "\n";
    return logpost;
}

//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
  gsl_vector_set_zero(z_hat);
  cHyp.pve=0.0;
  cHyp.pge=1.0;		
  return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	
	gsl_blas_ddot (Xb, Xb, &d);
	cHyp.pve=d/(double)ni_test;
	cHyp.pve/=cHyp.pve+1.0;
	cHyp.pge=1.0;
	
	gsl_vector_memcpy (z_hat, Xb);

	return;
}


void BSLMM::MCMC_Free_WorkVar(gsl_matrix *Result_hyp, gsl_matrix *Result_gamma, gsl_vector *z_hat)
{
    gsl_matrix_free(Result_hyp);
	gsl_matrix_free(Result_gamma);
	gsl_vector_free(z_hat);
    return;
}

void BSLMM::CreateSnpPosVec(vector<snpPos> &snp_pos)
{
    size_t pos;
    string rs;
    string chr;
    long int bp;
    size_t tt=0;
    
    for (size_t i=0; i < ns_total; ++i){
        if(!indicator_snp[i]) {continue;}

        rs = snpInfo[i].rs_number;
        chr = snpInfo[i].chr;
        bp = snpInfo[i].base_position;
        pos = tt;
        snpPos snp_temp={pos, rs, chr, bp};
        snp_pos.push_back(snp_temp);
        tt++;
    }
    
    snpInfo.clear(); // clear snpInfo vector
    
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // sort snp by chr and bp
    
    return;
}

void BSLMM::InitialMap(const vector<pair<size_t, double> > &pos_loglr, const vector<snpPos> &snp_pos){
    
    //with sorted vectors as input
    size_t pos;
    
    for (size_t i=0; i<ns_test; ++i) {
		mapRank2pos[i]=pos_loglr[i].first;
        mapPos2Rank[pos_loglr[i].first] = i;
        
        mapOrder2pos[i] = snp_pos[i].pos;
        mapPos2Order[snp_pos[i].pos] = i;
	}
    
    for (size_t i=0; i<ns_test; ++i) {
        pos = mapRank2pos[i];
		mapRank2Order[i] = mapPos2Order[pos];
        
        pos = mapOrder2pos[i];
        mapOrder2Rank[i] = mapPos2Rank[pos];
	}
    return;
}

void BSLMM::CreateGammaProposal(const double *p_gamma){
	//calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
	gsl_rng_env_setup();
	const gsl_rng_type * gslType;
	gslType = gsl_rng_default;
    
	if (randseed<0)
	{
		time_t rawtime;
		time (&rawtime);
		tm * ptm = gmtime (&rawtime);
		randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
	}
	gsl_r = gsl_rng_alloc(gslType);
	gsl_rng_set(gsl_r, randseed);

	gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma); // set up proposal function for gamma
}


double BSLMM::MHPropose(uchar **X_Genotype, const double *p_gamma, gsl_vector *z, const double &ztz, LModel &model_old, LModel &model_new, int &flag_gamma, double &logPost_new, double &logPost_old)
{
    int repeat=1;
    
    if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
    else {repeat=1;}
    
    double logMHratio=0.0;
    logMHratio+=ProposeHnRho(model_old.cHyp, model_new.cHyp, repeat);
    logMHratio+=ProposePi(model_old.cHyp, model_new.cHyp, repeat);
   // cout << "logMHratio from proposing pi = " << logMHratio << "\n";
    cout << "propose h, rho, pi success, proposing gamma..." << endl;
    
    logMHratio+=ProposeGamma (model_old, model_new, p_gamma, repeat, X_Genotype, z, ztz, flag_gamma); //JY
    cout << "logMHratio from proposing h, pi, gamma = " << logMHratio << "\n";
    
    if(flag_gamma==1) {
        nadd++;
        //cout << "add a snp\n";
    }
    else if(flag_gamma==2) {
        ndel++;
        //cout << "delete a snp\n";
    }
    else if(flag_gamma==3) {
        nswitch++;
       // cout << "switch a snp\n";
    }
    else nother++;
    //cout << "old rank = " ;
    //PrintVector(model_old.rank);
    //cout << "\n propose gamma success... with rank_new = " ;
   // PrintVector(model_new.rank);
    
    model_new.AssignVar(X_Genotype, z, mapRank2pos, ns_test);
   // cout << "Xtz_sub vector:";
  //  PrintVector(&model_new.Xty_sub.vector);
    
    if (model_new.rank.size() > 0 && model_new.rank.size() <= 5)
        { model_new.AssignVar(X_Genotype, z, mapRank2pos, ns_test); }
    else { SetXgamma (model_old, model_new, X_Genotype, z);
       // cout << "XtX_sub set from model_old:\n";
       // PrintMatrix(model_new.XtX, 3, model_new.rank.size());
       // model_new.AssignVar(X_Genotype, z, mapRank2pos, ns_test);
        //cout << "XtX_sub set from AssignVar:\n";
        //PrintMatrix(model_new.XtX, 3, model_new.rank.size());
    }
    logPost_new = model_new.CalcPosterior(ztz, ni_test, ns_test, gsl_r, a_mode, trace_G);

   // cout << "Calcposterior success." << endl;
    
    cout << "logPost_old = " << logPost_old << "; logPost_new = " << logPost_new << "\n";
    logMHratio+=logPost_new-logPost_old;
    cout << "final logMHratio = " << logMHratio << endl;

    return(logMHratio);

}


void BSLMM::MHmove(const bool &accept, const int &flag_gamma, LModel &model_old, LModel &model_new, double &logPost_old, double &logPost_new)
{
    if (accept==1) {
        
        if(flag_gamma==1) nadd_accept++;
        else if(flag_gamma==2) ndel_accept++;
        else if(flag_gamma==3) nswitch_accept++;
        else nother_accept++;
        
      // cout << "accept proposal ... " << endl;
    
        logPost_old=logPost_new;
        model_old.Copy(model_new);
    }
    else {
        model_new.cHyp = model_old.cHyp;
       // cout << "reject proposal ... " << endl;

    }
    
    return;
}

void BSLMM::MHsave(const size_t &t, size_t &w, gsl_matrix *Result_hyp, gsl_matrix *Result_gamma, LModel &model_old, vector<pair<double, double> > &beta_g, const double &mean_z, const vector<pair<size_t, double> > &pos_loglr, const vector<snpPos> &snp_pos, double &logPost_old)
{
    size_t w_col, pos;
    
    if (t%r_pace==0) {
        w_col=w%w_pace;
        if (w_col==0) {
            if (w==0) {
               // WriteResult (0, Result_hyp, Result_gamma, w_col);
            }
            else {
               // WriteResult (1, Result_hyp, Result_gamma, w_col);
                gsl_matrix_set_zero (Result_hyp);
                gsl_matrix_set_zero (Result_gamma);
            }
        }
        
        gsl_matrix_set (Result_hyp, w_col, 0, model_old.cHyp.h);
        gsl_matrix_set (Result_hyp, w_col, 1, model_old.cHyp.pve);
        gsl_matrix_set (Result_hyp, w_col, 2, model_old.cHyp.rho);
        gsl_matrix_set (Result_hyp, w_col, 3, model_old.cHyp.pge);
        gsl_matrix_set (Result_hyp, w_col, 4, model_old.cHyp.logp);
        gsl_matrix_set (Result_hyp, w_col, 5, model_old.cHyp.n_gamma);
        gsl_matrix_set (Result_hyp, w_col, 6, logPost_old);
        
        for (size_t i=0; i<model_old.cHyp.n_gamma; ++i) {
            
            pos=mapRank2pos[model_old.rank[i]]+1;
            gsl_matrix_set (Result_gamma, w_col, i, pos); //snp position starting from 1
            beta_g[pos-1].first+=gsl_vector_get(model_old.beta, i);
            beta_g[pos-1].second+=1.0;	
        }
        
        if (a_mode==13) {
            pheno_mean+=mean_z;
        }
        w++;
    }
    return;
}

//if a_mode==13, then run probit model
void BSLMM::MCMC (uchar **X_Genotype, gsl_vector *z) {
    
	clock_t time_start;
    
    double time_set=0, time_post=0; r_pace=1;
    
    LModel model_old, model_new;
    model_old.InitialVar(ni_test, s_max);
    model_new.InitialVar(ni_test, s_max);
    cout << "create two LModel structures success...\n";
    
    gsl_matrix *Result_hyp=gsl_matrix_alloc (w_pace, 7);
	gsl_matrix *Result_gamma=gsl_matrix_alloc (w_pace, s_max);
	gsl_vector *z_hat=gsl_vector_alloc (ni_test);

	//for quantitative traits, z is centered already in gemma.cpp, but just in case
	double mean_z = CenterVector (z);
    double ztz;
	gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    cout << "ztz = " << ztz << endl;
    cout << "mean of z" << mean_z << endl;

	double logPost_new, logPost_old, logMHratio;
	
	if (a_mode==13) {
		pheno_mean=0.0;
	}
	
	vector<pair<double, double> > beta_g;
	for (size_t i=0; i<ns_test; i++) {
		beta_g.push_back(make_pair(0.0, 0.0));
	}
	
	vector<pair<size_t, double> > pos_loglr;
	time_start=clock();
    cout << "start calculating marginal LRT...\n";
	//MatrixCalcLmLR (X_Genotype, z, pos_loglr, ns_test, ni_test, trace_G); //Simple linear regression save time?
    cout << "trace_G = trace(X'X) = " << trace_G;
    
    stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio
	time_Proposal=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

    // Jingjing add a vector of "snpPos" structs snp_pos
    vector<snpPos> snp_pos;
    CreateSnpPosVec(snp_pos);
   // printSNPInfo(snp_pos, 10); // print out a few markers' info
	
    InitialMap(pos_loglr, snp_pos);
    double *p_gamma = new double[ns_test];
	CalcPgamma (p_gamma); // calculate discrete distribution for gamma
	CreateGammaProposal(p_gamma);
	
	//initial parameters
    cout << "start initializing mcmc...\n";
    
	InitialMCMC (X_Genotype, z, model_old, pos_loglr, snp_pos);
	cHyp_initial=model_old.cHyp;
    model_old.AssignVar(X_Genotype, z, mapRank2pos, ns_test);
    cout << "print out first submatrix of XtX:\n";
    PrintMatrix(model_old.Xgamma, 10, 3);
    PrintMatrix(model_old.XtX, 3, model_old.rank.size());
    PrintVector(model_old.Xty);
    PrintVector(model_old.rank);
    size_t markeri_pos = mapRank2pos[1];
    snp_pos[markeri_pos].printMarker();
    
    
    model_old.getSubVar(model_old.rank.size());
    logPost_old = model_old.CalcPosterior(ztz, ni_test, ns_test, gsl_r, a_mode, trace_G);
    cout << "Initial logPost = " << logPost_old << endl;

	//calculate centered z_hat, and pve
	if (a_mode==13) {
		if (model_old.cHyp.n_gamma==0) {
			CalcCC_PVEnZ (z_hat, model_old.cHyp);
		}
		else {
			CalcCC_PVEnZ (model_old.Xbeta, z_hat, model_old.cHyp);
		}
	}
	
	//start MCMC
	bool accept;
	size_t total_step=w_step+s_step;
	size_t w=0, w_col; // pos declared earlier (JY)
    int flag_gamma=0;
    double accept_percent;
	
    cout << "start MCMC... \n" ;
	for (size_t t=0; t<total_step; ++t) {
        
		if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC (with acceptance percentage)", t, total_step-1, (double)n_accept/(double)(t*n_mh+1)); cout << "\n";}
//		if (t>10) {break;}
        
		if (a_mode==13) {			
			SampleZ (z, z_hat, z);
			mean_z=CenterVector (z);
			gsl_blas_ddot(z,z,&ztz);
					
            //First proposal
            model_old.AssignVar(X_Genotype, z, mapRank2pos, ns_test);
            model_old.getSubVar(model_old.rank.size());
            logPost_old = model_old.CalcPosterior(ztz, ni_test, ns_test, gsl_r, a_mode, trace_G);
		}

		//MH steps
		for (size_t i=0; i<n_mh; ++i) {
            
            time_start=clock();
            
            //cout << "start MHPropose... \n" ;
            logMHratio = MHPropose(X_Genotype, p_gamma, z, ztz, model_old, model_new, flag_gamma, logPost_new, logPost_old);
            
            time_set+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
            
            if (logMHratio>=0.0 || log(gsl_rng_uniform(gsl_r))<=logMHratio) {accept=1; n_accept++;}
			else {accept=0;}
            
            //cout << "start MHMove... \n" ;
            MHmove(accept, flag_gamma, model_old, model_new, logPost_old, logPost_new);
            
        }
        
		//calculate z_hat, and pve
		if (a_mode==13) {
			if (model_old.cHyp.n_gamma==0) {
				CalcCC_PVEnZ (z_hat, model_old.cHyp);
			}
			else {
				CalcCC_PVEnZ (model_old.Xbeta, z_hat, model_old.cHyp);
			}
			//sample mu and update z hat
			gsl_vector_sub (z, z_hat);
			mean_z+=CenterVector(z);
			mean_z+=gsl_ran_gaussian(gsl_r, sqrt(1.0/(double) ni_test) );			
			gsl_vector_add_constant (z_hat, mean_z);
		}
        
        
		
		//Save data
		if (t<w_step) {continue;}
		else {
            
            accept_percent = (double)n_accept/(double)(t*n_mh);
            if (accept_percent<0.000001) {
                cerr << "acceptance percentage = " << accept_percent << " < 0.000001; ABORT MCMC...";
                exit(1);
            }
            
           // cout << "start MHSave... \n" ;
			MHsave(t, w, Result_hyp, Result_gamma, model_old, beta_g, mean_z, pos_loglr, snp_pos, logPost_old);
		}
        
	}
        
	cout<<endl;
	cout<<"time on selecting Xgamma: "<<time_set<<endl;
	cout<<"time on calculating posterior: "<<time_post<<endl;

	w_col=w%w_pace;
	// WriteResult (1, Result_hyp, Result_gamma, w_col);
	
	gsl_vector *alpha=gsl_vector_alloc (ns_test);
	gsl_vector_set_zero (alpha);
	WriteParam (beta_g, alpha, w, snp_pos, pos_loglr);
	gsl_vector_free(alpha);
    
    model_old.FreeMem();
    model_new.FreeMem();
    MCMC_Free_WorkVar(Result_hyp, Result_gamma, z_hat);
    
	delete [] p_gamma;
	beta_g.clear();
	
	return;
}
