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
    n_type = cPar.n_type;
    mFunc = cPar.mFunc;
    e = cPar.e;
    vscale = cPar.vscale;
    FIXHYP = cPar.FIXHYP;
    saveSNP = cPar.saveSNP;
    iniType = cPar.iniType;
    iniSNPfile = cPar.iniSNPfile;
    hypfile = cPar.hypfile;
    
    UnCompBufferSize = cPar.UnCompBufferSize;
    CompBuffSizeVec = cPar.CompBuffSizeVec;
    Compress_Flag = cPar.Compress_Flag;
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

bool comp_pi (pair<string, double> a, pair<string, double> b)
{
    return (a.second > b.second);
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

void BSLMM::WriteIniSNP (const vector< pair<string, double> > &pivec, size_t n_snp)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".iniSNP";
    cout << "write iniSNP at " << file_str << endl;
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for (size_t i=0; i<n_snp; ++i) {
        outfile<< pivec[i].first <<endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
}

void BSLMM::WriteIniSNP (const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".iniSNP";
    //cout << "write iniSNP at " << file_str << endl;
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    for (size_t i=0; i<rank.size(); ++i) {
        outfile<< snp_pos[SNPrank_vec[rank[i]].second].rs <<endl;
    }
    
    outfile.clear();
    outfile.close();
    return;
}

void BSLMM::WriteMatrix(const gsl_matrix * X, const string filename){
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

void BSLMM::WriteVector(const gsl_vector * X, const string filename){
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


void PrintVector(const gsl_vector * x){
    for(size_t i=0; i < x->size; ++i){
        cout<<setprecision(6) << gsl_vector_get(x, i) << ", ";
    }
    cout << endl; 
}

void PrintVector(const gsl_vector * x, const size_t s){
    for(size_t i=0; i < s; ++i){
        cout<<setprecision(6) << gsl_vector_get(x, i) << ", ";
    }
    cout << endl;
}

void PrintMatrix(const gsl_matrix * X, const size_t nrow, const size_t ncol){
    for (size_t i=0; i<nrow; i++) {
        gsl_vector_const_view row = gsl_matrix_const_subrow(X, i, 0, ncol);
        PrintVector(&row.vector);
    }
}

void PrintVector(const vector <double> &x){
    for(size_t i=0; i<x.size(); ++i){
        cout <<setprecision(6) << x[i] << ", ";
    }
    cout << endl; 
}

void PrintVector(const vector <double> &x, const size_t s){
    for(size_t i=0; i<s; ++i){
        cout <<setprecision(6) << x[i] << ", ";
    }
    cout << endl;
}

void PrintVector(const vector <size_t> &x){
    for(size_t i=0; i<x.size(); ++i){
        cout <<setprecision(6) << x[i] << ", ";
    }
    cout << endl;
}

void PrintVector(const double *x){
    for(size_t i=0; i<41; ++i){
        cout <<setprecision(6) << x[i] << ", ";
    }
    cout << endl; 
}

void PrintVector(const  uchar *x, const size_t length){
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

void BSLMM::WriteParamtemp(vector<pair<double, double> > &beta_g, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr)
{
    string file_str;
    file_str="./output/"+file_out;
    file_str+=".paramtemp";
    
    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
    
    //outfile<<"markerID"<<"\t"<<"chr"<<"\t" <<"bp"<<"\t" << "maf" << "\t" << "Func_code"<< "\t" <<"beta"<<"\t"<<"gamma" << endl;
    
    size_t pos;
    vector< pair<string , double> > pivec;
    string rs;
    double pitemp;
    
    for (size_t i=0; i<ns_test; ++i) {
        
        // save the data along the order of all variants, snp_pos is sorted by order
        rs = snp_pos[i].rs;
        outfile<< rs <<"\t"<< snp_pos[i].chr<<"\t" <<snp_pos[i].bp << "\t";
        
        for (size_t j=0; j < n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                outfile << j << "\t";
                break;
            }
            else if(j == (n_type - 1)) outfile << "NA" << "\t";
        }
        
        outfile << scientific << setprecision(6)  << snp_pos[i].maf << "\t";
        
        pos = snp_pos[i].pos;
        if (beta_g[pos].second!=0) {
            pitemp = beta_g[pos].second/(double)s_step;
            outfile << beta_g[pos].first/beta_g[pos].second<< "\t" << pitemp <<endl;
        }
        else {
            pitemp = 0.0;
            outfile << 0.0 << "\t" << 0.0 << endl;
        }
        pivec.push_back(make_pair(rs, pitemp));
    }
    
    // Save top significant SNPs as starting position for next step
    /*if (saveSNP) {
        size_t n_snp = 0;
        stable_sort(pivec.begin(), pivec.end(), comp_pi);
        double q_genome = gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
        
        for (size_t i=0; i<pos_loglr.size(); ++i) {
            if (2.0*pos_loglr[i].second>q_genome) {n_snp++;}
        }

        if (n_snp<10) {n_snp = 10;}
        if (n_snp>s_max) {n_snp = s_max;}
        WriteIniSNP(pivec, n_snp);
    }*/
    
    outfile.clear();	
    outfile.close();
}


void BSLMM::WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr)
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"geno_pos"<<"\t"<<"markerID"<<"\t" << "maf" <<"\t"<<"chr"<<"\t"
			<<"bp"<<"\t" << "Func_code"  << "\t"<<"alpha"<<"\t" << "lrt" << "\t"
			<<"beta"<<"\t"<<"gamma" << endl; //JY added gamma_var
	
    
    size_t r, pos;
	for (size_t i=0; i<ns_test; ++i) {
        
        // save the data along the order of all variants, snp_pos is sorted by order
		outfile<<snp_pos[i].pos << "\t" << snp_pos[i].rs<<"\t" << snp_pos[i].maf  <<"\t"<< snp_pos[i].chr<<"\t"
		<<snp_pos[i].bp<<"\t";
        
        for (size_t j=0; j < n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                outfile << j << "\t";
                continue;
            }
            else if (j == (n_type - 1) ) {
                outfile << "NA" << "\t";
            }
        }
		
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
        else {cerr << "pos dose not match snp_pos[i].pos...\n"; exit(-1);}
        
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

void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col)
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
        
        outfile_hyp<<"pve \t pge \t n_gamma \t logPosterior \t h \t rho0 \t rho1 \t theta0 \t theta1 \t sigma0 \t sigma 1"<<endl;
        
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
            outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, 0)<<"\t" <<gsl_matrix_get (Result_hyp, i, 1) << "\t" ;
            outfile_hyp << (int)gsl_matrix_get (Result_hyp, i, 2)<< "\t" <<scientific << setprecision(6) << gsl_matrix_get (Result_hyp, i, 3)<< "\t"<< gsl_matrix_get (Result_hyp, i, 4)<< "\t"<< gsl_matrix_get (Result_hyp, i, 5)<< "\t"<< gsl_matrix_get (Result_hyp, i, 6)<< "\t"<< gsl_matrix_get (Result_hyp, i, 7)<< "\t"<< gsl_matrix_get (Result_hyp, i, 8)<< "\t"<< gsl_matrix_get (Result_hyp, i, 9) << "\t"<< gsl_matrix_get (Result_hyp, i, 10)<<endl;
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
			outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, 0)<<"\t" <<gsl_matrix_get (Result_hyp, i, 1) << "\t";
            outfile_hyp << (int)gsl_matrix_get (Result_hyp, i, 2)<< "\t" <<scientific << setprecision(6) << gsl_matrix_get (Result_hyp, i, 3)<<endl;
		}
		
		for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<s_max; ++j) {
				outfile_gamma<<(int)gsl_matrix_get (Result_gamma, i, j)<<"\t";
			}
			outfile_gamma<<endl;
		}
        
        for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<n_type; ++j) {
				outfile_theta<<scientific<<setprecision(6) << exp(gsl_matrix_get (Result_theta, i, j))<<"\t";
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
		pos=SNPrank_vec[rank[i]].first;
		gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
        getGTgslVec(X, &Xgamma_col.vector, pos, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
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

//InitialMCMC with LModel input
/*void BSLMM::InitialMCMC (uchar **X, const gsl_vector *Uty, LModel &model, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos)

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

 gsl_vector_free(Xtxvec);
 gsl_matrix_free(Xr);
 gsl_matrix_free(XtXr);
 gsl_vector_free(Xtyr);
 gsl_vector_free(xvec);
 gsl_vector_free(yres);
     
     
     stable_sort (model.rank.begin(), model.rank.end(), comp_vec); //sort the initial rank.
     PrintVector(model.rank);
     
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
 } */


void BSLMM::setHyp(double htemp, double theta_temp, double subvar_temp){
    
    h = htemp;
    
    theta.assign(n_type, theta_temp);
    log_theta.clear();
    log_theta.push_back(log(theta[0]));
    log_theta.push_back(log(theta[1]));
    
    subvar.assign(n_type, subvar_temp);
    
    
    if (hypfile.c_str() == NULL) {
        cout << "no hypefile ... ";
        exit(1);
    }
    else{
        ifstream infile(hypfile.c_str(), ifstream::in);
        if(!infile) {cout << "Error opening file " << hypfile << endl; exit(-1);
           // cout << "load hyp from hypfile... " << hypfile << endl;
    }
    string line;
    char *pch, *nch;
    //cout << "load fixed hyper parameter values from : " << hypfile << endl;
    
    while (!safeGetline(infile, line).eof()) {
        if ((line[0] == 'h') || (line[0] == '#') || (line[0] == 'p')) {
            continue;
        }
        else{
            pch = (char *)line.c_str();
            nch = strchr(pch, '\t');
            h = strtod(pch, NULL);
            pch = (nch == NULL) ? NULL : nch+1;
            //cout << "h = "<<setprecision(5) << h ;
            
            nch = strchr(pch, '\t');
            theta[0] = strtod(pch, NULL);
            log_theta[0] = log(theta[0]);
            pch = (nch == NULL) ? NULL : nch+1;
            
            nch = strchr(pch, '\t');
            theta[1] = strtod(pch, NULL);
            log_theta[1] = log(theta[1]);
            pch = (nch == NULL) ? NULL : nch+1;
            
            nch = strchr(pch, '\t');
            subvar[0] = tau * strtod(pch, NULL);
            pch = (nch == NULL) ? NULL : nch+1;
            
            subvar[1] = tau * strtod(pch, NULL);
        }
    }
    }
    cout << "pve = "<<setprecision(5) << h << "; theta = " << theta[0] << ", " << theta[1] << "; subvar = " << subvar[0] << ", " << subvar[1] << endl;
}


//InitialMCMC currently used
void BSLMM::InitialMCMC (uchar **X, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos)

{
    double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
    cHyp.n_gamma=0;
    for (size_t i=0; i<pos_loglr.size(); ++i) {
        if (2.0*pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
    }
    //cout << "number of snps before adjust = " << cHyp.n_gamma << endl;
    if (cHyp.n_gamma<10) {cHyp.n_gamma=10;}
    if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
    if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}
    
    
    if (!iniSNPfile.empty() && iniType == 0) {
        
        ifstream infile(iniSNPfile.c_str(), ifstream::in);
        if(!infile) {cout << "Error opening file " << iniSNPfile << endl; exit(-1);}
        string lineid;
        rank.clear();
        size_t orderj, rankj;
        
        cout << "Start loading initial snp IDs from " << iniSNPfile << "\n";
        
        while (!safeGetline(infile, lineid).eof()) {
            
            orderj = 0;
            for (size_t i=0; i < snp_pos.size(); i++) {
                if (snp_pos[i].rs.compare(lineid) == 0) {
                    orderj=i;
                    rankj = SNPorder_vec[orderj].second;
                    rank.push_back(rankj);
                    //cout << lineid << " with rank = " << rankj;
                    //snp_pos[orderj].printMarker();
                    break;
                }
                if ( i == (snp_pos.size()-1) ) {
                    cout << "Reach end of vector, failed finding snp " << lineid << endl;
                }
            }
            
        }
        cout << endl;
        infile.close();
        infile.clear();
        cHyp.n_gamma = rank.size();
    }
    else if(iniType == 0) {iniType = 1;}
    else if(iniType == 1) {
        cout << "Start with top SVT SNPs.\n";
        rank.clear();
        for (size_t i=0; i<cHyp.n_gamma; ++i) {
            rank.push_back(i);
        }
    } // Start with most significant variants from SVT
    else if(iniType == 2) {
        cout << "Start with top SVT SNPs/ColinearTest.\n";
        size_t posr, j=0, i=0;
        double xtx;
        
        rank.clear();
        rank.push_back(0);
        posr = SNPrank_vec[0].first;
        xtx = XtX_diagvec[posr];
        // cout << "rank added: " << 0 << ", ";
        
        gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
        gsl_vector * xvec = gsl_vector_alloc(ni_test);
        getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); //get geno column
        
        gsl_matrix * XtXr = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
        gsl_matrix * XtXr_temp = gsl_matrix_alloc(cHyp.n_gamma, cHyp.n_gamma);
        gsl_vector * Xtxvec = gsl_vector_alloc(cHyp.n_gamma);
        
        do{
       // for (size_t i=1; i < cHyp.n_gamma; ++i){
            gsl_matrix_set_col(Xr, i, xvec);
            gsl_matrix_set(XtXr, i, i, xtx);
            
            if (i>0) {
                gsl_matrix_const_view Xr_sub = gsl_matrix_const_submatrix(Xr, 0, 0, ni_test, i);
                gsl_vector_view Xtxvec_sub = gsl_vector_subvector(Xtxvec, 0, i);
                gsl_blas_dgemv(CblasTrans, 1.0, &Xr_sub.matrix, xvec, 0.0, &Xtxvec_sub.vector);
                
                gsl_vector_view XtX_subrow = gsl_matrix_subrow(XtXr, i, 0, i);
                gsl_vector_view XtX_subcol = gsl_matrix_subcolumn(XtXr, i, 0, i);
                gsl_vector_memcpy(&XtX_subrow.vector, &Xtxvec_sub.vector);
                gsl_vector_memcpy(&XtX_subcol.vector, &Xtxvec_sub.vector);
                
            }
            gsl_matrix_memcpy(XtXr_temp, XtXr);
            
            if ((rank.size() < cHyp.n_gamma)) {
            
                do{
                    j++; // Consider rank j
                    //cout << "consider rank j" << j << endl;
                } while( (j < ns_test) && ColinearTest(X, Xr, XtXr_temp, j, rank.size()));
                
                rank.push_back(j);
                posr = SNPrank_vec[j].first;
                xtx = XtX_diagvec[posr];
                getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); // get geno column
            }
            i++;
            
        }while(i < (cHyp.n_gamma));
        
        //PrintMatrix(XtXr, cHyp.n_gamma, cHyp.n_gamma);
        
        gsl_matrix_free(Xr);
        gsl_matrix_free(XtXr);
        gsl_matrix_free(XtXr_temp);
        gsl_vector_free(xvec);
        gsl_vector_free(Xtxvec);

    } // Start with most significant variants from SVT
    else if(iniType == 3){
        cout << "Start with Step-wise SNPs. \n";
    vector<pair<size_t, double> > rank_loglr;
    size_t posr, radd;
        
    for (size_t i=1; i<1000; ++i) {
        rank_loglr.push_back(make_pair(i, pos_loglr[i].second));
    }
    cout << endl;
    
    rank.clear();
    rank.push_back(0);
    posr = SNPrank_vec[0].first;
   // cout << "rank added: " << 0 << ", ";

    gsl_matrix * Xr = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
    gsl_matrix_set_zero(Xr);
    gsl_vector * xvec = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); //get geno column
    
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
            posr = SNPrank_vec[rank_loglr[j].first].first;
            getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); // get geno column
            rank_loglr[j].second = CalcLR(yres, xvec, posr);
        }
        stable_sort (rank_loglr.begin(), rank_loglr.end(), comp_lr); //sort the initial rank.
        
        radd = rank_loglr[0].first;
       // cout << "rank added: " << radd << ", ";
        posr = SNPrank_vec[radd].first;
        getGTgslVec(X, xvec, posr, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
        rank.push_back(radd);
        rank_loglr.erase(rank_loglr.begin());
    }
    
    gsl_matrix_free(Xr);
    gsl_matrix_free(XtXr);
    gsl_vector_free(Xtyr);
    gsl_vector_free(xvec);
    gsl_vector_free(yres);
    gsl_vector_free(Xtxvec);
    }
    cout << "number of snps = " << cHyp.n_gamma << endl;
    stable_sort (rank.begin(), rank.end(), comp_vec); //sort the initial rank.
    PrintVector(rank);
    
    cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
    cHyp.h=pve_null;
    
    if (cHyp.logp==0) {cHyp.logp=-0.000001;}
    if (cHyp.h==0) {cHyp.h=0.1;}
    
    gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
    SetXgamma (UtXgamma, X, rank);
   // cout<<"initial value of h = "<< h <<endl;
    
    double sigma_a2;
    if (trace_G!=0) {
        sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp));
    } else {
        sigma_a2=cHyp.h*1.0/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
    }
    if (sigma_a2==0) {sigma_a2=0.025;}
   // cout << "initial sigma_a2 = " << sigma_a2 << endl;
    
    //cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
    gsl_matrix_free (UtXgamma);
    
    if (cHyp.rho>1.0) {cHyp.rho=1.0;}
    if (cHyp.h<h_min) {cHyp.h=h_min;}
    if (cHyp.h>h_max) {cHyp.h=h_max;}
    //if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
    //if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
    if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
    if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}
    
    //cout << "start setHyp... \n";
    setHyp(cHyp.h, ((double)cHyp.n_gamma/(double)ns_test), sigma_a2);
    cHyp.theta = theta;
    cHyp.log_theta = log_theta;
    cHyp.subvar = subvar; // initial subvar vector
    
   // e_shape =  e;
   // e_rate = e_shape / sigma_a2; // Gamma with mean sigma_a2,
    //cout << "IG with shape = " << e_shape << "; rate = " << e_rate;
   // cout << "; mean = " << sigma_a2 << "; sd = " << (sigma_a2 / sqrt(e_shape)) << endl;
    
    
    cout<<"initial value of h = "<< h <<endl;
    //cout<<"initial value of rho = "<<cHyp.rho<<endl;
    cout<<"initial value of theta_vec = "<<theta[0] << ", " << theta[1] <<endl;
    cout << "initial value of sub-variance_vec = " << subvar[0]<< ", " << subvar[1] << endl;
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
/*double BSLMM::CalcPosterior (const double yty, class HYPBSLMM &cHyp, const gsl_vector *pi_vec, const vector<size_t> &rank)
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
    // calculate gamma likelihood
    logpost += CalcLikegamma(pi_vec, rank);
    
    return logpost;
}*/

// for the new model
double BSLMM::CalcPosterior (const double yty, class HYPBSLMM &cHyp)
{
    double logpost=0.0;
    
    //for quantitative traits, calculate pve and pge
    //pve and pge for case/control data are calculted in CalcCC_PVEnZ
    if (a_mode==11) {
        cHyp.pve=0.0;
     //   cHyp.pge=1.0;
    }
    //calculate likelihood
    if (a_mode==11) {logpost-=0.5*(double)ni_test*log(yty);}
    else {logpost-=0.5*yty;}
    // calculate gamma likelihood
    //logpost += CalcLikegamma(cHyp);
    
    return logpost;
}

// exponentiate log_theta vector
void expVector(vector<double> &expvec, vector<double> &logvec){
    expvec.clear();
    for (size_t i=0; i < logvec.size(); i++) {
        expvec.push_back(exp(logvec[i]));
    }
}

void BSLMM::CalcPivec(const vector<double> &theta, gsl_vector *pi_vec, const vector<SNPPOS> &snp_pos){
    double pi_temp, p = 1.0 / double(ns_test);
    gsl_vector_set_all(pi_vec, p);
   // cout << "theta vec: "; PrintVector(theta);
    // Pivec in the same order as snp_pos
    for (size_t i=0; i < ns_test; i++) {
        pi_temp = 0.0;
       // if (i < 300) {
         //   cout << "weight: ";
          //  PrintVector(snp_pos[i].weight);
        //}
        
        for (size_t j=0; j < n_type; j++) {
            
            if(snp_pos[i].indicator_func[j])
            {
                if ( (snp_pos[i].weight[j] <= 0.0) || (snp_pos[i].weight[j] > 1.0) ) {
                    cerr << "ERROR: weight[j] = " << snp_pos[i].weight[j] << endl;
                    exit(-1);
                }
                else{
                    pi_temp += theta[j] * snp_pos[i].weight[j];
                }
            }
        }
        if (pi_temp >1 || pi_temp < 0) {
            cerr << "pi_temp = " << pi_temp << endl;
            exit(-1);
        }
        else gsl_vector_set(pi_vec, i, pi_temp);
       // if (i < 300) {
         //   cout << "pi" << i << " = " << pi_temp << endl;
       // }
    }
}

void BSLMM::CalcSvec(const vector<double> &subvar, gsl_vector *sigma_vec, const vector<SNPPOS> &snp_pos){
    
    double sigma_temp;
    gsl_vector_set_zero(sigma_vec);
    for (size_t i=0; i < ns_test; i++) {
        sigma_temp=0.0;
        for (size_t j=0; j < n_type; j++)
        {
            if(snp_pos[i].indicator_func[j])
                sigma_temp += subvar[j] * snp_pos[i].weight[j];
        }
        if (sigma_temp <= 0 ) {
            cerr << "Error: sigma"<< i<< " = " << sigma_temp << "<= 0" << endl;
            exit(-1);
        }
        else gsl_vector_set(sigma_vec, i, sigma_temp);
    }
}

void BSLMM::getSubVec(gsl_vector *sigma_subvec, const gsl_vector * sigma_vec, const vector<size_t> &rank)
{
    size_t order_i;
    for (size_t i=0; i < rank.size(); i++) {
        order_i = mapRank2Order[rank[i]];
        gsl_vector_set(sigma_subvec, i, gsl_vector_get(sigma_vec, order_i));
    }
}

//calculate subvar vector given log_theta, h and rho
void BSLMM::setSubvar(class HYPBSLMM &cHyp, const vector<double> &Gvec)
{
    for (size_t i=0; i < n_type; i++) {
        
        cHyp.subvar[i] = cHyp.h * cHyp.rho_vec[i] / (1.0 - cHyp.h);
        
       // cHyp.subvar[i] = (double)ni_test * cHyp.h * cHyp.rho_vec[i] / ((1.0 - cHyp.h) * Gvec[i] * theta[i]);
        
        //cHyp.subvar[i] = (double)ni_test * cHyp.h * cHyp.rho_vec[i] / ((1.0 - cHyp.h) * Gvec[i] );
       
       /* if (cHyp.m_gamma[i] > 1) {
            cHyp.subvar[i] = (double)mFunc[i] * cHyp.h * cHyp.rho_vec[i] / ((1.0 - cHyp.h) * (double)cHyp.m_gamma[i] * Gvec[i] );
        }
        else cHyp.subvar[i] = (double)mFunc[i] * cHyp.h * cHyp.rho_vec[i] / ((1.0 - cHyp.h) * Gvec[i] ); */
    }
}

void BSLMM::setSubvar(class HYPBSLMM &cHyp)
{
    for (size_t i=0; i < n_type; i++) {
        cHyp.subvar[i] = cHyp.h * cHyp.rho_vec[i] / (1.0 - cHyp.h);
    }
}

//set sigma_subvec vectoer.
void BSLMM::getSubVec(gsl_vector *sigma_subvec, class HYPBSLMM &cHyp, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    size_t order_i;
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = mapRank2Order[rank[i]];
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                gsl_vector_set(sigma_subvec, i, cHyp.subvar[j]);
                continue;
            }
        }
    }
}

//Used for EM-Block
void BSLMM::getSubVec(gsl_vector *sigma_subvec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    size_t order_i;
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = SNPrank_vec[rank[i]].second;
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                gsl_vector_set(sigma_subvec, i, subvar[j]);
                continue;
            }
        }
    }
}

//set sigma_subvec and mgamma vectoer and trace vector Gvec

void BSLMM::set_mgamma(class HYPBSLMM &cHyp, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    size_t order_i;
    
    cHyp.m_gamma.assign(n_type, 0);
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = SNPrank_vec[rank[i]].second;
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                cHyp.m_gamma[j]++;
                continue;
            }
        }
    }
}

void BSLMM::setGvec(class HYPBSLMM &cHyp, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos, const gsl_matrix * XtX, vector<double> &Gvec)
{
    size_t order_i;

    cHyp.m_gamma.assign(n_type, 0);
    Gvec.assign(n_type, 0.0);
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = mapRank2Order[rank[i]];
        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[order_i].indicator_func[j]) {
                cHyp.m_gamma[j]++;
                Gvec[j] += gsl_matrix_get(XtX, i, i);
                continue;
            }
        }
    }
}


void BSLMM::CalcSvec(const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos){
    
    if (cHyp.n_gamma != rank.size()) {
        cerr << "Error: cHyp.n_gamma not equal to the size of rank\n";
        exit(-1);
    }
    size_t order_i;
    double sigma_temp;
    gsl_vector_set_zero(sigma_vec);
    for (size_t i=0; i < cHyp.n_gamma; i++) {
        order_i = mapRank2Order[rank[i]];
        sigma_temp = CalcSigma(cHyp, order_i, snp_pos);
        gsl_vector_set(sigma_vec, i, sigma_temp);
    }
}

double BSLMM::CalcSigma(const class HYPBSLMM &cHyp, const size_t &order_i, const vector<SNPPOS> &snp_pos){
    double sigma_temp = 0.0;
    for (size_t j=0; j < n_type; j++)
    {
        if(snp_pos[order_i].indicator_func[j])
                sigma_temp += cHyp.subvar[j] * snp_pos[order_i].weight[j];
    }
    if (sigma_temp <= 0 ) {
        cerr << "Error: sigma = " << sigma_temp << "<= 0" << endl;
        exit(-1);
    }
    return sigma_temp;
}

// Calculate prior P(subvar_j)
double BSLMM::CalcPsubvar (const class HYPBSLMM &cHyp, size_t j)
{
    double logp_var;
    logp_var = (e_shape - 1.0) * log(cHyp.subvar[j]) - e_rate * (cHyp.subvar[j]);
    return logp_var;
}

// Calculate prior P(theta_j)
double BSLMM::CalcPtheta (const class HYPBSLMM &cHyp, size_t j)
{
    return (-cHyp.log_theta[j]);
}

// Calculate prior P(1/subvar) , gamma(shape, rate)
double BSLMM::CalcPsubvar (const class HYPBSLMM &cHyp)
{
    double logp_var = 0.0, sum_logv = 0.0, sum_v = 0.0;
    for (size_t i=0; i<n_type; i++) {
        sum_logv += log(cHyp.subvar[i]);
        sum_v += e_rate * cHyp.subvar[i];
    }
    logp_var = (e_shape - 1.0) * sum_logv - sum_v ;
    return (logp_var);
}

// Calculate prior P(theta)
double BSLMM::CalcPtheta (const class HYPBSLMM &cHyp)
{
    double logp_theta = 0.0;
    for (size_t i=0; i<n_type; i++) {
        logp_theta += cHyp.log_theta[i];
    }
    return (-logp_theta);
}


// Calculate P(gamma | theta)
double BSLMM::CalcLikegamma(const gsl_vector *pi_vec, const vector<size_t> &rank)
{
    double pi_rank=0.0, pi_nonrank=0.0;
    size_t order_i, s_size = rank.size();
    
    for (size_t i=0; i < ns_test; i++) {
        pi_nonrank += log(1.0 - gsl_vector_get(pi_vec, i));
    }
    for (size_t i=0; i < s_size; i++) {
        order_i = mapRank2Order[rank[i]];
        pi_rank += log(gsl_vector_get(pi_vec, order_i));
        pi_nonrank -= log(1.0 - gsl_vector_get(pi_vec, order_i));
    }
    return (pi_rank + pi_nonrank);
}

/*double BSLMM::CalcLikegamma(const class HYPBSLMM &cHyp)
{
    double loglikegamma = 0.0;
    
    for (size_t i=0; i < n_type; i++) {
        loglikegamma +=  cHyp.log_theta[i] * (double)cHyp.m_gamma[i] + (double)(mFunc[i] - cHyp.m_gamma[i]) * log(1.0 - exp(cHyp.log_theta[i]));
    }
    return loglikegamma;
}*/

//used in EM_BLock
double BSLMM::CalcLikegamma(const class HYPBSLMM &cHyp)
{
    double loglikegamma = 0.0;
    
    for (size_t i=0; i < n_type; i++) {
        loglikegamma +=  log_theta[i] * (double)cHyp.m_gamma[i] + (double)(mFunc[i] - cHyp.m_gamma[i]) * log(1.0 - theta[i]);
    }
    return loglikegamma;
}


void CalcXVbeta(gsl_matrix *X, const gsl_vector * sigma_vec)
{
    for (size_t i=0; i < X->size1; i++) {
        gsl_vector_view xvec = gsl_matrix_row(X, i);
        gsl_vector_scale(&xvec.vector, gsl_vector_get(sigma_vec, i));
    }
}


/*double BSLMM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *pi_vec, gsl_vector *sigma_vec, const vector<size_t> &rank)
{
    //conditioning on hyper parameters: subvar, log_theta
    
    double logpost=0.0;
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = rank.size();
    
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);
    gsl_vector *Xty_temp=gsl_vector_alloc (s_size);
        
    //calculate Omega
    gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
    CalcXVbeta(Omega, &sigma_sub.vector);
    gsl_matrix_set_identity (M_temp);
    gsl_matrix_add (Omega, M_temp);
    
    //calculate beta_hat
    gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
    logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);	//solve Omega * beta_hat = Xty for beta_hat
    // Omega was inverted here
    // logdet_0 = det(Omega)
    //cout << "inv(Omega) * Xty: "; PrintVector(beta_hat);
    gsl_vector_mul(beta_hat, &sigma_sub.vector);
    //cout << "beta_hat: "; PrintVector(beta_hat);
    
    gsl_blas_ddot (Xty_temp, beta_hat, &d);
    P_yy-=d;
    if (P_yy <= 0) {
        cerr << "Error in calcposterior: P_yy = " << P_yy << endl;
        
        cout << "h = " << cHyp.h << "; rho: " << cHyp.rho_vec[0] << ", " << cHyp.rho_vec[1];
        cout << "; theta: " << exp(cHyp.log_theta[0]) << ", " << exp(cHyp.log_theta[1]);
        cout << "; subvar: " << cHyp.subvar[0] << ", " << cHyp.subvar[1];
        cout << "beta_hat: "; PrintVector(beta_hat);
        
        cerr << "set beta_hat to 0\n";
        gsl_vector_set_zero(beta_hat);
        P_yy = yty;
        //
        //exit(-1);
    }
    
    //sample tau
    double tau=1.0;
    if (a_mode==11) {tau = gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/P_yy); }
    
    //sample beta
    for (size_t i=0; i<s_size; i++)
    {
        gsl_vector_set(beta, i, gsl_ran_gaussian(gsl_r, 1));
    }
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, &beta_sub.vector);
    double beta_i, sigma_i;
    for (size_t i=0; i<s_size; i++) {
        beta_i = gsl_vector_get(&beta_sub.vector, i);
        sigma_i = gsl_vector_get(&sigma_sub.vector, i);
        //cout << "beta_i = " << beta_i << "; sigma_i = " << sigma_i << "; tau = " << tau<< endl;
        gsl_vector_set(&beta_sub.vector, i, beta_i * sqrt(sigma_i/tau));
    }
    //cout << "Sampled beta with mean 0: "; PrintVector(&beta_sub.vector);
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
    //cout << "-0.5 * logdet_0 = " << logpost;
    
    if (a_mode==11) {logpost-=0.5*(double)ni_test*log(P_yy);}
    else {logpost-=0.5*P_yy;}
    
   // cout << "; log_P(Y | ..) = " << logpost ;
    double loglikegamma = CalcLikegamma(pi_vec, rank);
    //Calc gamma likelihood P(gamma | theta) for posterior
    logpost += loglikegamma;
   // cout << " + loglikegamma = " << loglikegamma << "; logpost = " << logpost << endl;

    gsl_matrix_free (Omega);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (Xty_temp);

    return logpost;
} */

// for the new model
double BSLMM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag)
{
    //conditioning on hyper parameters: subvar, log_theta
    
    double logpost=0.0;
    double d;
    double P_yy=0.0, logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;
    
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
    gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
    gsl_matrix *Omega_temp=gsl_matrix_alloc (s_size, s_size);
    gsl_vector *beta_hat=gsl_vector_alloc (s_size);
    
    //calculate Omega
    gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
   // cout << "Omega : "; PrintMatrix(Omega, 5, 5);
    CalcXVbeta(Omega, &sigma_sub.vector);
    gsl_vector_view Omega_diag = gsl_matrix_diagonal(Omega);
    gsl_vector_add_constant(&Omega_diag.vector, 1);
    gsl_matrix_memcpy(Omega_temp, Omega);
    
    //calculate beta_hat
    //cout << "inv(Omega) * Xty: ";
    double lambda = 0.0;
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(Omega, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.00001;
    
    
    /*gsl_permutation * p = gsl_permutation_alloc(s_size);
    int s = s_size;
    gsl_linalg_LU_decomp(Omega, p, &s);
    gsl_linalg_LU_solve(Omega, p, &Xty_sub.vector, beta_hat);
    gsl_permutation_free(p);*/
    
    //logdet_O = LapackCholSolve(Omega_temp, &Xty_sub.vector, beta_hat);
    //cout << "lambda = " << lambda << endl;
    //EigenSolve(Omega, &Xty_sub.vector, beta_hat, lambda);
    //EigenSolve(Omega, &Xty_sub.vector, beta_hat);
    if(LapackSolve(Omega, &Xty_sub.vector, beta_hat)!=0)
       EigenSolve(Omega, &Xty_sub.vector, beta_hat);
    logdet_O=LapackLogDet(Omega_temp);
    //logdet_O = CalcLogdet(Omega_temp);
    
    //cout << "beta_hat from solve : "; PrintVector(beta_hat);
    gsl_vector_mul(beta_hat, &sigma_sub.vector);
    //cout << "beta_hat: "; PrintVector(beta_hat);
    
    double bxy;
    gsl_blas_ddot (&Xty_sub.vector, beta_hat, &bxy);
    
    if (yty < bxy) {
        
        do{
            if (s_size < 10) {
                // Error_Flag = 1;
                //cerr << "Error in calcPosterior: P_yy = " << P_yy << endl;
                //cout << "P_yy = " << P_yy << endl;
                //cout << "est beta_hat: "; PrintVector(beta_hat);
                WriteMatrix(&Xgamma_sub.matrix, "_X");
                WriteMatrix(&XtX_sub.matrix, "_XtX");
                WriteMatrix(Omega, "_Omega");
                WriteVector(&sigma_sub.vector, "_sigma");
                WriteVector(&Xty_sub.vector, "_Xty");
                WriteVector(beta_hat, "_beta");
            }

         gsl_vector_add_constant(&Omega_diag.vector, lambda);
         if(LapackSolve(Omega, &Xty_sub.vector, beta_hat)!=0)
             EigenSolve(Omega, &Xty_sub.vector, beta_hat);
         gsl_blas_ddot (&Xty_sub.vector, beta_hat, &bxy);
        } while (yty < bxy) ;
    }
    
    gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
    gsl_vector_memcpy(&beta_sub.vector, beta_hat);
    gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);
    gsl_blas_ddot (Xb, Xb, &d);
    if (a_mode==11) {
        cHyp.pve=d/(double)ni_test;
        //cHyp.pve/=cHyp.pve+1.0/tau;
        // cHyp.pge=1.0;
    }
    
    logpost = -((double)cHyp.n_gamma * logrv + logdet_O) + tau * bxy;
    logpost *= 0.5;
    
    P_yy = (yty - bxy);
    
    
    gsl_matrix_free (Omega_temp);
    gsl_matrix_free (Omega);
    gsl_vector_free (beta_hat);
    
    return logpost;
}

// for the new model
//calculate likelihood P(Y | gamma, subvar, theta)
double BSLMM::CalcLikelihood (const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const class HYPBSLMM &cHyp, gsl_vector *sigma_vec, bool &Error_Flag)
{
    //double sigma_a2=cHyp.h/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test * trace_G);
    
    double loglike=0.0;
    double d, P_yy=yty, logdet_O=0.0;
    size_t s_size = cHyp.n_gamma;
    
  if (s_size == 0) {
        //calculate likelihood if ngamma=0
        if (a_mode==11) {loglike-=0.5*(double)ni_test*log(yty);}
        else {loglike-=0.5*yty;}
  }
    
  else{
   // gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);
    gsl_vector_const_view sigma_sub = gsl_vector_const_subvector(sigma_vec, 0, s_size);
    
      gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
      gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
      gsl_vector *beta_hat=gsl_vector_alloc (s_size);
      gsl_vector *Xty_temp=gsl_vector_alloc (s_size);
      
      //calculate Omega
      gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
      CalcXVbeta(Omega, &sigma_sub.vector);
      gsl_matrix_set_identity (M_temp);
      gsl_matrix_add (Omega, M_temp);
      
      //calculate beta_hat
      gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);
      logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);	//solve Omega * beta_hat = Xty for beta_hat
      // Omega was inverted here
      // logdet_0 = det(Omega)
      gsl_vector_mul(beta_hat, &sigma_sub.vector);
      gsl_blas_ddot (Xty_temp, beta_hat, &d);
      P_yy-=d;
      if (P_yy <= 0) {
          Error_Flag = 1;
          cout << "Error in calclikelihood: P_yy = " << P_yy << endl;
          cout << "h = "<<setprecision(6) << cHyp.h << "; rho: " << cHyp.rho_vec[0] << ", " << cHyp.rho_vec[1];
          cout << "theta: "<<setprecision(6) << exp(cHyp.log_theta[0]) << ", " << exp(cHyp.log_theta[1]);
          cout << "; subvar: "<<setprecision(6) << cHyp.subvar[0] << ", " << cHyp.subvar[1];
          cout << "beta_hat: "; PrintVector(beta_hat);
          
          cout << "set beta_hat to 0\n";
          gsl_vector_set_zero(beta_hat);
          P_yy = yty;
         // exit(-1);
      }
      else {Error_Flag = 0;}
    
    loglike = -0.5 * logdet_O;
    if (a_mode==11) {loglike -= 0.5 * (double)ni_test * log(P_yy);}
    else {loglike -= 0.5*P_yy;}
    
    gsl_matrix_free (Omega);
    gsl_matrix_free (M_temp);
    gsl_vector_free (beta_hat);
    gsl_vector_free (Xty_temp);
  }
    
    return loglike;
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
            if (rank_new.size() > 0) {
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0) || (ColinearTest(X, Xgamma_old, XtX_old, r_add, rank_new.size())));
            }
            else {
                do {
                    r_add=gsl_ran_discrete (gsl_r, gsl_t);
                } while ((mapRank2in.count(r_add)!=0));
            }
            
            
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
            o_remove = SNPrank_vec[r_remove].second;
            rank_new.erase(rank_new.begin()+col_id);
            size_t s_size = rank_new.size();
            mapRank2in.erase(r_remove);
            
            if (s_size > 0) {
              if (cHyp_new.n_gamma<=2 || cHyp_old.n_gamma<=2) {
                SetXgamma (Xgamma_temp, X, rank_new);
                CalcXtX (Xgamma_temp, z, s_size, XtX_gamma, Xtz_gamma);
              } else {
                SetXgamma (X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_temp, XtX_gamma, Xtz_gamma);
              }
                CalcRes(Xgamma_temp, z, XtX_gamma, Xtz_gamma, z_res, s_size, ztz);
                gsl_s = MakeProposal(o_remove, p_BFr, X, z_res, mapRank2in);
            }
            else {
                gsl_s = MakeProposal(o_remove, p_BFr, X, z, mapRank2in);
            }
            //construct gsl_s, JY
            do {
                j_add = gsl_ran_discrete(gsl_r, gsl_s);
                o_add = (o_remove - win) + j_add;
                r_add = SNPorder_vec[o_add].second;
            } while ((s_size > 0) && ColinearTest(X, Xgamma_temp, XtX_gamma, r_add, s_size));
            
            //cout << "o_add = " << o_add <<  "; r_add = "<<r_add << endl;
            if((o_add < 0) || (o_add >= (long int)ns_test) || (o_add == (long int)o_remove))
                cout << "ERROR proposing switch snp"; //new snp != removed snp
            
            gsl_a = MakeProposal(o_add, p_BFa, X, z_res, mapRank2in);
            
            double prob_total_remove=1.0;
            double prob_total_add=1.0;
            
            for (size_t ii=0; ii<rank_new.size(); ++ii) {
                r = rank_new[ii];
                o = SNPrank_vec[r].second;
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
           // cout << "successfully switched a snp" << endl;
        }
        
        else {logp+=0.0;}//do not change
        
    }
    
    stable_sort (rank_new.begin(), rank_new.end(), comp_vec);
    mapRank2in.clear();
    return logp;
}

void BSLMM::WriteHyptemp(gsl_vector *LnPost, vector<double> &em_gamma){
    
    //vector<double> em_gamma_avg(n_type, 0.0);
    //em_gamma_avg[0] = em_gamma[0] / (double)s_step;
    //em_gamma_avg[1] = em_gamma[1] / (double)s_step;
    
    double em_logpost = 0.0, logpost_max =  gsl_vector_max(LnPost);
    cout << "logpost_max = " << logpost_max << endl;
    for (size_t i=0; i < s_step; i++) {
        em_logpost += exp(gsl_vector_get(LnPost, i) - logpost_max);
    }
    em_logpost /= double(s_step);
    em_logpost = log(em_logpost) + logpost_max;
    
    //save E(GV, sumbeta2[0], sumbeta2[1], sum_m0, sum_m1, lnpost, m0, m1, n0, n1)
    string file_hyp;
    file_hyp = "./output/" + file_out;
    file_hyp += ".hyptemp";
    ofstream outfile_hyp;
    outfile_hyp.open (file_hyp.c_str(), ofstream::out);
    if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
    
    outfile_hyp << scientific << setprecision(6) << (GV / (double)s_step) << "\t" << rv << "\t" << sumbeta2[0] << "\t" << sumbeta2[1] << "\t";
    outfile_hyp << em_gamma[0] << "\t" << em_gamma[1]<< "\t";
    outfile_hyp << em_logpost << "\t" << (em_gamma[0] / (double)s_step) << "\t" << (em_gamma[1] / (double)s_step) << "\t" ;
    outfile_hyp << mFunc[0] << "\t"  << mFunc[1] << endl;
    outfile_hyp.clear();
    outfile_hyp.close();
    
}



void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col, const vector<SNPPOS> &snp_pos, const vector<pair<size_t, double> > &pos_loglr)
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


//Testing current version of MCMC

void BSLMM::MCMC_Test (uchar **X, const gsl_vector *y, bool original_method) {
    
    if (original_method) {
        cout << "Run Block MCMC...\n";
    }
    clock_t time_start;
    double time_set=0, time_post=0;
    cout << "# of unique function types = " << n_type << endl;
    
    //new model related
    //gsl_vector *pi_vec = gsl_vector_alloc (ns_test);
   // gsl_vector *sigma_vec = gsl_vector_alloc(ns_test);
   // gsl_vector_set_zero(sigma_vec);
    gsl_vector *sigma_subvec_old = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_old);
    gsl_vector *sigma_subvec_new = gsl_vector_alloc(s_max);
    gsl_vector_set_zero(sigma_subvec_new);
    gsl_vector *LnPost = gsl_vector_alloc(s_step); //save logPost...
    vector<double> em_gamma(n_type, 0.0); //save sum of m0, m1
    GV = 0.0; sumbeta2.assign(n_type, 0.0);
    
    //same as old model
    //gsl_matrix *Result_hyp=gsl_matrix_alloc (w_pace, 4);
    //gsl_matrix *Result_gamma=gsl_matrix_alloc (w_pace, s_max);
    
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
    double mean_z = CenterVector (z);
    gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    gsl_vector_scale(z, 1.0 / sqrt(ztz / (double)(ni_test-1))); // standardize phenotype z
    //for quantitative traits, y is centered already in gemma.cpp, but just in case
    gsl_blas_ddot(z, z, &ztz); // ztz is the sum square of total SST
    
    rv = 0.6; tau = 1.0 / rv;
    logrv = log(rv);
    cout << "ztz = " << ztz << "; Fix Residual Variance = " << rv << endl;
    cout << "PI = " << M_PI << "; tau = " << tau << "; logrv = " <<logrv << endl;
    
    // cout << "mean of z = " << mean_z << endl;
    
    //Initialize variables for MH
    double logPost_new, logPost_old;
    double logMHratio;
    vector<size_t> rank_new, rank_old;
    vector<double> Gvec;
    class HYPBSLMM cHyp_old, cHyp_new;
    bool Error_Flag=0;
    
    if (a_mode==13) {
        pheno_mean=0.0;
    }
    vector<pair<double, double> > beta_g; //save beta estimates
    for (size_t i=0; i<ns_test; i++) {
        beta_g.push_back(make_pair(0.0, 0.0));
    }
    
    //Setup log-likelihood ratio test statistics
    time_start=clock();
    
    // Jingjing add a vector of "snpPos" structs snp_pos
    vector<SNPPOS> snp_pos;
    CreateSnpPosVec(snp_pos); //ordered by chr/bp
    cout << "1 / SNP_sd  = "; PrintVector(SNPsd, 10);
    
    vector<pair<size_t, double> > pos_loglr;
    MatrixCalcLmLR (X, z, pos_loglr, ns_test, ni_test, SNPsd, Gvec, XtX_diagvec, snp_pos, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); //calculate trace_G or Gvec
    trace_G = (Gvec[0] + Gvec[1]) / double(ni_test * ns_test);
    cout << "trace_G = " << trace_G << endl;
    cout << "Total trace_G vec : " << Gvec[0] << ", " << Gvec[1] << endl;
    
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp);
    stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio
    
    //Jingjing's edit, create maps between rank and order
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
    
    SNPorder_vec.clear();
    SNPrank_vec.clear();
    for (size_t i=0; i<ns_test; i++) {
        SNPorder_vec.push_back(make_pair(snp_pos[i].pos, mapOrder2Rank[i]));
        SNPrank_vec.push_back(make_pair(pos_loglr[i].first, mapRank2Order[i]));
    }
    //end of Jingjing's edit
    
    //Calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t
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
    
    //Initial parameters
    cout << "Start initializing MCMC ... \n";
    InitialMCMC (X, z, rank_old, cHyp_old, pos_loglr, snp_pos); // Initialize rank and cHyp
    inv_subvar.assign(n_type, 0.0), log_subvar.assign(n_type, 0.0);
    inv_subvar[0] = (1.0 / subvar[0]); inv_subvar[1] = (1.0 / subvar[1]);
    log_subvar[0] = (log(subvar[0])); log_subvar[1] = (log(subvar[1]));
    cout << "inv_subvar = "; PrintVector(inv_subvar);
    cout << "log_subvar = "; PrintVector(log_subvar);
    
    if (cHyp_old.n_gamma > 0) {
        SetXgamma (Xgamma_old, X, rank_old);
        CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
    }
    
    //cout << "Set m_gamma... \n";
    set_mgamma(cHyp_old, rank_old, snp_pos);
    cout << "initial m_gamma: " << cHyp_old.m_gamma[0] << ", "<< cHyp_old.m_gamma[1]<< endl;
    cout << "Set sigma_subvec... \n";
    getSubVec(sigma_subvec_old, rank_old, snp_pos);
    //PrintVector(sigma_subvec_old, rank_old.size());
    
    cHyp_initial=cHyp_old;
    gsl_vector_memcpy(sigma_subvec_new, sigma_subvec_old);
    
    //Calculate first loglikelihood
    //cout << "first calculating logpost ... \n";
    if (cHyp_old.n_gamma==0) {
        logPost_old = CalcPosterior (ztz, cHyp_old) + CalcLikegamma(cHyp_old);
    }
    else {
        logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag) + CalcLikegamma(cHyp_old);
    }
    if (!Error_Flag) {
       // cout <<  "logPost_old = " << logPost_old << endl;
    }
    else {
        cerr << "Failed at initialMCMC...\n";
        exit(-1);
    }
    cout <<  "Initial logPost_old = " << logPost_old << endl;
    
    /*
    cout << "First 10 responses: \n";
    PrintVector(z, 10);
    cout << "First 10 Xtz: \n";
    PrintVector(Xtz_old, 10);
    cout << "First 10 X : \n";
    PrintMatrix(Xgamma_old, 10, 10);
    cout << "First 10 XtX : \n";
    PrintMatrix(XtX_old, 10, 10);
     */
    
    //calculate centered z_hat, and pve
    if (a_mode==13) {
        if (cHyp_old.n_gamma==0) {
            CalcCC_PVEnZ (z_hat, cHyp_old);
        }
        else {
            CalcCC_PVEnZ (Xb_old, z_hat, cHyp_old);
        }
    }
    
    //Start MCMC
    w_pace=1000;
    int accept; // accept_theta; naccept_theta=0,
    size_t total_step=w_step+s_step;
    size_t repeat=1;
    int flag_gamma=0;
    double accept_percent, betai; // accept_theta_percent;
    
    cHyp_new = cHyp_old;
    rank_new = rank_old;

    for (size_t t=0; t<total_step; ++t) {
        
       // if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));}
        //		if (t>10) {break;}
        
        if (a_mode==13) {
            SampleZ (y, z_hat, z); //sample z
            mean_z=CenterVector (z);
            gsl_blas_ddot(z,z,&ztz);
            
            //First proposal need to be revised
            if (cHyp_old.n_gamma==0) {
                logPost_old= CalcPosterior (ztz, cHyp_old) + CalcLikegamma(cHyp_old);
            } else {
                gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
                gsl_vector_view Xtz_sub=gsl_vector_subvector(Xtz_old, 0, rank_old.size());
                gsl_blas_dgemv (CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0, &Xtz_sub.vector);
                logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag) + CalcLikegamma(cHyp_old);
            }
        }

        //////// Set repeat number
        //if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
        //else {repeat=1;}
        //cout << "n_mh = " << n_mh << endl;
        
        for (size_t i=0; i<n_mh; ++i) {

           // cout << "propose gamam...\n";
           // cout << "old rank: "; PrintVector(rank_old);
            //repeat = 1;
            logMHratio = ProposeGamma (rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat, X, z, Xgamma_old, XtX_old, Xtz_old,  ztz, flag_gamma); //JY
            //rank_new.clear(); cHyp_new.n_gamma=0;
            if(flag_gamma==1) nadd++;
            else if(flag_gamma==2) ndel++;
            else if(flag_gamma==3) nswitch++;
            else nother++;
           // cout << "new rank: "; PrintVector(rank_new);
            //cout << "flag_gamma = " << flag_gamma << endl;
            //cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
           // cout << "propose gamma MHratio = " << exp(logMHratio) << endl;
            
            if (rank_new.size() > 0) {
                //this if makes sure that rank_old.size()==rank_remove.size() does not happen
               // cout << "start set Xgamma_new... " << endl;
                if (cHyp_new.n_gamma<=2 || cHyp_old.n_gamma<=2) {
                    SetXgamma (Xgamma_new, X, rank_new);
                    CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
                } else {
                    SetXgamma (X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_new, XtX_new, Xtz_new);
                }
                set_mgamma(cHyp_new, rank_new, snp_pos);
                getSubVec(sigma_subvec_new, rank_new, snp_pos);
            }
            else {
                cHyp_new.m_gamma.assign(n_type, 0);
            }
           // cout << "new m_gamma: " << cHyp_new.m_gamma[0] << ", "<< cHyp_new.m_gamma[1]<< endl;
            if (rank_new.size()==0) {
                logPost_new = CalcPosterior (ztz, cHyp_new) + CalcLikegamma(cHyp_new);
            } else {
                logPost_new = CalcPosterior (Xgamma_new, XtX_new, Xtz_new, ztz, Xb_new, beta_new, cHyp_new, sigma_subvec_new, Error_Flag) + CalcLikegamma(cHyp_new);
            }
              //cout << "Calcposterior success." << endl;
            if (!Error_Flag) {
                logMHratio += logPost_new-logPost_old;
                // cout <<"logPost_old = " << logPost_old<< "; logPost_new = "<< logPost_new<< "; MHratio = " << exp(logMHratio) << endl;
                if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio)
                    { accept=1; if (flag_gamma < 4) n_accept++;}
                else {accept=0;}
            }
            else{
                accept=0;
            }
            
            //cout << "accept = " << accept << endl;
            //accept = 1;
            
            if (accept==1) {
                    if(flag_gamma==1) nadd_accept++;
                    else if(flag_gamma==2) ndel_accept++;
                    else if(flag_gamma==3) nswitch_accept++;
                    else nother_accept++;
                
                    logPost_old=logPost_new;
                    cHyp_old.n_gamma = cHyp_new.n_gamma;
                    cHyp_old.m_gamma = cHyp_new.m_gamma;
                    cHyp_old.pve = cHyp_new.pve;
                    //cHyp_old.rv = cHyp_new.rv;
                    //cHyp_old.pge = cHyp_new.pge;
                    gsl_vector_memcpy (Xb_old, Xb_new);
                    rank_old = rank_new;
                if(rank_new.size()>0){
                    gsl_vector_view sigma_oldsub=gsl_vector_subvector(sigma_subvec_old, 0, rank_old.size());
                    gsl_vector_view sigma_newsub=gsl_vector_subvector(sigma_subvec_new, 0, rank_old.size());
                    gsl_vector_memcpy(&sigma_oldsub.vector, &sigma_newsub.vector);
                
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
                cHyp_new.n_gamma = cHyp_old.n_gamma;
                cHyp_new.m_gamma = cHyp_old.m_gamma;
                rank_new = rank_old;
            }
            //cout << "copy data from new propose -> old " << endl;
        } //end of n_mh
        
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
        
         //if (t % 10 == 0 && t > w_step) {
         if (t % w_pace == 0 && t > w_step) {
             accept_percent = (double)n_accept/(double)((t+1) * n_mh);
             cout << "cHyp_old.n_gamma= " << cHyp_old.n_gamma << endl;
            cout << "gamma acceptance percentage = " <<setprecision(6) << accept_percent << endl ;
             cout << "m_gamma: " << cHyp_old.m_gamma[0] << ", " << cHyp_old.m_gamma[1]<< endl;
             cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); cout << endl;
             cout << "loglike: " << logPost_old << endl;
        }
        
        //Save data
        if (t<w_step) {continue;}
        else {
                gsl_vector_set (LnPost, (t-w_step), logPost_old);
                em_gamma[0] += (double)cHyp_old.m_gamma[0];
                em_gamma[1] += (double)cHyp_old.m_gamma[1];
                GV += cHyp_old.pve;
        
            if (cHyp_old.n_gamma > 0){
                for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
                    // beta_g saved by position
                    pos=SNPrank_vec[rank_old[i]].first;
                    betai = gsl_vector_get(beta_old, i);
                    beta_g[pos].first += betai;
                    beta_g[pos].second += 1.0;
                    for (size_t j=0; j < n_type; j++) {
                        if (snp_pos[SNPrank_vec[rank_old[i]].second].indicator_func[j]) {
                            sumbeta2[j] += betai * betai;
                            break;
                        }
                    }
                }
              }
            else{cout << "rank size = " << rank_old.size() << endl;}
            }
    }
    
    cout<< "MCMC completed ... " << endl << endl;
    accept_percent = (double)n_accept/(double)(total_step * n_mh);
    cout << "gamma acceptance percentage = " << accept_percent << endl ;
    cout << "m_gamma: " << cHyp_old.m_gamma[0] << ", " << cHyp_old.m_gamma[1]<< endl;
    cout << "beta_hat: "; PrintVector(beta_old, rank_old.size()); cout << endl;
    cout << "loglike: " << logPost_old << endl;
    
    //save last causal SNPIDs
    if (saveSNP) WriteIniSNP(rank_old, snp_pos);
    
    //Save temp EM results
    WriteHyptemp(LnPost, em_gamma);
    WriteParamtemp(beta_g, snp_pos, pos_loglr);
    
   // gsl_matrix_free(Result_hyp);
   // gsl_matrix_free(Result_gamma);
    
    gsl_vector_free(sigma_subvec_old);
    gsl_vector_free(sigma_subvec_new);
    gsl_vector_free(LnPost);
    
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
    
    delete [] p_gamma;
    beta_g.clear();
    
    return;
}
// end of current version









//acceptance percentage 15% MCMC version
/*
 void BSLMM::MCMC (uchar **X, const gsl_vector *y, bool original_method) {
 
 if (original_method) {
 cout << "Run previous working version of MCMC...\n";
 }
 clock_t time_start;
 double time_set=0, time_post=0;
 cout << "# of unique function types = " << n_type << endl;
 
 //new model related
 //gsl_vector *pi_vec = gsl_vector_alloc (ns_test);
 // gsl_vector *sigma_vec = gsl_vector_alloc(ns_test);
 // gsl_vector_set_zero(sigma_vec);
 gsl_vector *sigma_subvec_old = gsl_vector_alloc(s_max);
 gsl_vector_set_zero(sigma_subvec_old);
 gsl_vector *sigma_subvec_new = gsl_vector_alloc(s_max);
 gsl_vector_set_zero(sigma_subvec_new);
 
 //same as old model
 gsl_matrix *Result_hyp=gsl_matrix_alloc (w_pace, 11);
 gsl_matrix *Result_gamma=gsl_matrix_alloc (w_pace, s_max);
 //gsl_matrix_set_zero (Result_gamma);
 
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
 // cout << "ztz = " << ztz << endl;
 // cout << "mean of z = " << mean_z << endl;
 
 //Initialize variables for MH
 double logPost_new, logPost_old;
// double logPostTheta_new, logPostTheta_old;
// double logPostHRho_new, logPostHRho_old, loglikegamma, loglikebeta;
     double logMHratio;// logMHratio_theta, logMHratio_HRho;
 vector<size_t> rank_new, rank_old;
 vector<double> Gvec, Gvec_old, Gvec_new;
 class HYPBSLMM cHyp_old, cHyp_new;
 bool Error_Flag;
 
 if (a_mode==13) {
 pheno_mean=0.0;
 }
 vector<pair<double, double> > beta_g;
 for (size_t i=0; i<ns_test; i++) {
 beta_g.push_back(make_pair(0.0, 0.0));
 }
 
 //Setup log-likelihood ratio test statistics
 time_start=clock();
 
 // Jingjing add a vector of "snpPos" structs snp_pos
 vector<SNPPOS> snp_pos;
 CreateSnpPosVec(snp_pos); //ordered by chr/bp
 
 vector<pair<size_t, double> > pos_loglr;
 MatrixCalcLmLR (X, z, pos_loglr, ns_test, ni_test, Gvec, snp_pos, CompBuffSizeVec, UnCompBufferSize, Compress_Flag); //calculate trace_G or Gvec
 cout << "Total trace_G/n : " << Gvec[0] << ", " << Gvec[1] << endl;
 
 stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp);
 stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr); // sort log likelihood ratio
 
 //Jingjing's edit, create maps between rank and order
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
 iniType=1; // Start with top significant SNP by SVT
 
 // theta = {}; subvar = {}; //Initialize hyper parameters for theta and subvar
 InitialMCMC (X, z, rank_old, cHyp_old, pos_loglr, snp_pos); // Initialize rank and cHyp
 // subvar.clear(); subvar.push_back(); subvar.push_back();
 if (cHyp_old.n_gamma > 0) {
 SetXgamma (Xgamma_old, X, rank_old);
 CalcXtX (Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
 }
 
 cout << "Set trace vector and m_gamma... \n";
 setGvec(cHyp_old, rank_old, snp_pos, XtX_old, Gvec_old);
 cout << "Trace_G given gamma : " << Gvec_old[0] << ", " << Gvec_old[1] << endl;
 cout << "initial m_gamma: " << cHyp_old.m_gamma[0] << ", "<< cHyp_old.m_gamma[1]<< endl;
 
 cout << "Calculate sigma vectors... \n";
 setSubvar(cHyp_old, Gvec);
 cout << "initial subvar: " << cHyp_old.subvar[0] << ", "<< cHyp_old.subvar[1] << endl;
 getSubVec(sigma_subvec_old, rank_old, snp_pos);
 
 cHyp_initial=cHyp_old;
 gsl_vector_memcpy(sigma_subvec_new, sigma_subvec_old);
 
 //Calculate first loglikelihood
 if (cHyp_old.n_gamma==0) {
 logPost_old = CalcPosterior (ztz, cHyp_old) + CalcLikegamma(cHyp_old) + CalcPrho(cHyp_old);
 //logPostHRho_old = logPost_old + CalcPtheta(cHyp_old) + CalcPrho(cHyp_old);
 // loglikegamma = CalcLikegamma(cHyp_old);
 // logPostHRho_old = logPost_old - loglikegamma;
 // logPostTheta_old = CalcLikegamma(cHyp_old) + CalcPtheta(cHyp_old);
 }
 else {
 // loglikebeta = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old);
 logPost_old = CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag) + CalcLikegamma(cHyp_old) + CalcPrho(cHyp_old);
 // logPostHRho_old = logPost_old + CalcPtheta(cHyp_old) + CalcPrho(cHyp_old);
 // loglikegamma = CalcLikegamma(cHyp_old);
 // logPostHRho_old = logPost_old - loglikegamma;
 // logPostTheta_old = CalcLikegamma(cHyp_old) + CalcPtheta(cHyp_old);
 }
 if (!Error_Flag) {
 cout <<  "logPost_old = " << logPost_old << endl;
     //<<  "; logPostHRho_old = " << logPostHRho_old << endl;
 }
 else {
 cerr << "Failed at initialMCMC...\n";
 exit(-1);
 }*/
 
 // cout << "; logPostTheta_old = " << logPostTheta_old << endl;
 /*
 cout << "First 10 responses: \n";
 PrintVector(z, 10);
 cout << "First 10 Xtz: \n";
 PrintVector(Xtz_old, 10);
 cout << "First 10 X : \n";
 PrintMatrix(Xgamma_old, 10, 10);
 cout << "First 10 XtX : \n";
 PrintMatrix(XtX_old, 10, 10);
 */ 
/*
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

     int accept;// accept_hrho, accept_theta;
//size_t naccept_theta=0, naccept_hrho=0;
size_t total_step=w_step+s_step;
size_t w=0, w_col; // pos declared earlier (JY)
//	size_t repeat=0;
size_t repeat=1;
int flag_gamma=0;
     double accept_percent; // accept_hrho_percent, accept_theta_percent;

// double d_logp = min(0.5, (logp_max-logp_min)*logp_scale);
// cout << "d_log = " << d_logp << "; logp_min = " << logp_min << "; logp_max" << logp_max << endl;
cHyp_new = cHyp_old;
rank_new = rank_old;
cout << "logp_max = " << logp_max << "; logp_min = "<< logp_min ;
cout << "; logp_stepsize = " << (logp_max-logp_min)*logp_scale<< endl;

for (size_t t=0; t<total_step; ++t) {
    if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));}
    //		if (t>10) {break;}
    if (a_mode==13) {
        SampleZ (y, z_hat, z); //sample z
        mean_z=CenterVector (z);
        gsl_blas_ddot(z,z,&ztz);
        
        //First proposal need to be revised
        if (cHyp_old.n_gamma==0) {
            logPost_old=CalcPosterior (ztz, cHyp_old) + CalcLikegamma(cHyp_old) + CalcPrho(cHyp_old);
        } else {
            gsl_matrix_view Xold_sub=gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
            gsl_vector_view Xtz_sub=gsl_vector_subvector(Xtz_old, 0, rank_old.size());
            gsl_blas_dgemv (CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0, &Xtz_sub.vector);
            logPost_old=CalcPosterior (Xgamma_old, XtX_old, Xtz_old, ztz, Xb_old, beta_old, cHyp_old, sigma_subvec_old, Error_Flag) + CalcLikegamma(cHyp_old) + CalcPrho(cHyp_old);
        }
    }
    
    // set repeat number
    //if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
    //else {repeat=1;}
    
    //sample theta from beta distribution
    cHyp_old.log_theta[0] = log(gsl_ran_beta(gsl_r, 1 + cHyp_old.m_gamma[0], 3 + mFunc[0] - cHyp_old.m_gamma[0]));
    cHyp_old.log_theta[1] = log(gsl_ran_beta(gsl_r, 1 + cHyp_old.m_gamma[1], 3 + mFunc[1] - cHyp_old.m_gamma[1]));
    
    // sample h, rho, gamma
    for (size_t i=0; i<n_mh; ++i) {
        // cout << "propose gamam...\n";
        // cout << "old rank: "; PrintVector(rank_old);
        logMHratio = ProposeHnRho(cHyp_old, cHyp_new, repeat);
        logMHratio += ProposeGamma (rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat, X, z, Xgamma_old, XtX_old, Xtz_old,  ztz, flag_gamma); //JY
        stable_sort (rank_new.begin(), rank_new.end(), comp_vec);
        if(flag_gamma==1) nadd++;
        else if(flag_gamma==2) ndel++;
        else if(flag_gamma==3) nswitch++;
        else nother++;
        // cout << "new rank: "; PrintVector(rank_new);
        //cout << "flag_gamma = " << flag_gamma << endl;
        //cout << "propose gamma success... with rank_new.size = " << rank_new.size() << endl;
        // cout << "propose gamma MHratio = " << exp(logMHratio) << endl;
        
        if (cHyp_new.n_gamma > 0) {
            //this if makes sure that rank_old.size()==rank_remove.size() does not happen
            if (cHyp_new.n_gamma<=20 || cHyp_old.n_gamma<=20) {
                SetXgamma (Xgamma_new, X, rank_new);
                CalcXtX (Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
            } else {
                //cout << "start set Xgamma_new... " << endl;
                SetXgamma (X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new, Xgamma_new, XtX_new, Xtz_new);
            }
        }
        setGvec(cHyp_new, rank_new, snp_pos, XtX_new, Gvec_new);
        setSubvar(cHyp_new, Gvec);
        getSubVec(sigma_subvec_new, rank_new, snp_pos);
        // cout << "new m_gamma: " << cHyp_new.m_gamma[0] << ", "<< cHyp_new.m_gamma[1]<< endl;
        if (cHyp_new.n_gamma==0) {
            logPost_new=CalcPosterior (ztz, cHyp_new) + CalcLikegamma(cHyp_new) + CalcPrho(cHyp_new);
        } else {
            logPost_new=CalcPosterior (Xgamma_new, XtX_new, Xtz_new, ztz, Xb_new, beta_new, cHyp_new, sigma_subvec_new, Error_Flag) + CalcLikegamma(cHyp_new) + CalcPrho(cHyp_new);
        }
        //  cout << "Calcposterior success." << endl;
        if (!Error_Flag) {
            logMHratio += logPost_new-logPost_old;
            // cout <<"logPost_old = " << logPost_old<< "; logPost_new = "<< logPost_new<< "; MHratio = " << exp(logMHratio) << endl;
            if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio)
            { accept=1; if (flag_gamma < 4) n_accept++;}
            else {accept=0;}
        }
        else{
            accept=0;
        }
        
        // cout << "accept = " << accept << endl;
        
        if (accept==1) {
            if(flag_gamma==1) nadd_accept++;
            else if(flag_gamma==2) ndel_accept++;
            else if(flag_gamma==3) nswitch_accept++;
            else nother_accept++;
            
            logPost_old=logPost_new;
            Gvec_old = Gvec_new;
            cHyp_old.h = cHyp_new.h;
            cHyp_old.rho_vec = cHyp_new.rho_vec;
            cHyp_old.subvar = cHyp_new.subvar;
            cHyp_old.n_gamma = cHyp_new.n_gamma;
            cHyp_old.m_gamma = cHyp_new.m_gamma;
            cHyp_old.pve = cHyp_new.pve;
            cHyp_old.pge = cHyp_new.pge;
            gsl_vector_memcpy (Xb_old, Xb_new);
            rank_old = rank_new;
            gsl_vector_view sigma_oldsub=gsl_vector_subvector(sigma_subvec_old, 0, rank_old.size());
            gsl_vector_view sigma_newsub=gsl_vector_subvector(sigma_subvec_new, 0, rank_old.size());
            gsl_vector_memcpy(&sigma_oldsub.vector, &sigma_newsub.vector);
            
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
        } else {
            cHyp_new.h = cHyp_old.h;
            cHyp_new.rho_vec = cHyp_old.rho_vec;
            cHyp_new.subvar = cHyp_old.subvar;
            cHyp_new.n_gamma = cHyp_old.n_gamma;
            cHyp_new.m_gamma = cHyp_old.m_gamma;
            Gvec_new = Gvec_old;
            rank_new = rank_old;
        }
        //cout << "copy data from new propose -> old " << endl;
    } //end of n_mh
    
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
    
    
    // if (t % 10 == 0 && t > w_step) {
    if (t % w_pace == 0 && t > w_step) {
        // accept_theta_percent = (double)naccept_theta / (double)(t+1);
        //accept_hrho_percent = (double)naccept_hrho / (double)(t+1);
        accept_percent = (double)n_accept/(double)((t+1) * n_mh);
        
        // cout << "theta acceptance percentage = " << accept_theta_percent << endl ;
        //cout << "HRho acceptance percentage = " << accept_hrho_percent << endl ;
        cout << "gamma acceptance percentage = " << accept_percent << endl ;
        
        cout << "h = " << cHyp_old.h << "; rho: " << cHyp_old.rho_vec[0] << ", " << cHyp_old.rho_vec[1]<< endl;
        cout << "theta: " << exp(cHyp_old.log_theta[0]) << ", " << exp(cHyp_old.log_theta[1])<< endl;
        cout << "m_gamma: " << cHyp_old.m_gamma[0] << ", " << cHyp_old.m_gamma[1]<< endl;
        cout << "subvar: " << cHyp_old.subvar[0] << ", " << cHyp_old.subvar[1] << endl;
        cout << "sigma_subvec_old : " <<setprecision(6);
        PrintVector(sigma_subvec_old, 5);
        cout << "theta_est: "; PrintVector(theta);
        cout << "beta_hat: "; PrintVector(beta_old, 5); cout << endl;
    }
    
    //Save data
    if (t<w_step) {continue;}
    else {
        if (t%r_pace==0) {
            w_col=w%w_pace;
            if (w_col==0) {
                if (w==0) {
                    WriteResult (0, Result_hyp, Result_gamma, w_col);
                }
                else {
                    WriteResult (1, Result_hyp, Result_gamma, w_col);
                    gsl_matrix_set_zero (Result_hyp);
                    gsl_matrix_set_zero (Result_gamma);
                }
            }
            
            gsl_matrix_set (Result_hyp, w_col, 0, cHyp_old.pve);
            gsl_matrix_set (Result_hyp, w_col, 1, cHyp_old.pge);
            gsl_matrix_set (Result_hyp, w_col, 2, cHyp_old.n_gamma);
            gsl_matrix_set (Result_hyp, w_col, 3, logPost_old);
            gsl_matrix_set (Result_hyp, w_col, 4, cHyp_old.h);
            gsl_matrix_set (Result_hyp, w_col, 5, cHyp_old.rho_vec[0]);
            gsl_matrix_set (Result_hyp, w_col, 6, cHyp_old.rho_vec[1]);
            gsl_matrix_set (Result_hyp, w_col, 7, exp(cHyp_old.log_theta[0]));
            gsl_matrix_set (Result_hyp, w_col, 8, exp(cHyp_old.log_theta[1]));
            gsl_matrix_set (Result_hyp, w_col, 9, cHyp_old.subvar[0]);
            gsl_matrix_set (Result_hyp, w_col, 10, cHyp_old.subvar[1]);
            
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
cout<< " MCMC completed... " <<endl;
// accept_theta_percent = (double)naccept_theta / (double)(total_step+1);
//accept_hrho_percent = (double)naccept_hrho / (double)(total_step);
accept_percent = (double)n_accept/(double)(total_step * n_mh);
//cout << "theta acceptance percentage = " << accept_theta_percent << endl ;
//cout << "HRho acceptance percentage = " << accept_hrho_percent << endl ;
cout << "gamma acceptance percentage = " << accept_percent << endl ;
cout << "h = " << cHyp_old.h << "; rho: " << cHyp_old.rho_vec[0] << ", " << cHyp_old.rho_vec[1]<< endl;
cout << "theta: " << exp(cHyp_old.log_theta[0]) << ", " << exp(cHyp_old.log_theta[1])<< endl;
cout << "m_gamma: " << cHyp_old.m_gamma[0] << ", " << cHyp_old.m_gamma[1]<< endl;
cout << "subvar: " << cHyp_old.subvar[0] << ", " << cHyp_old.subvar[1] << endl;
cout << "sigma_subvec_old : " <<setprecision(6); PrintVector(sigma_subvec_old, 5);
cout << "beta_hat: "; PrintVector(beta_old, 5); cout << endl;

cout<<"time on selecting Xgamma: "<<time_set<<endl;
cout<<"time on calculating posterior: "<<time_post<<endl;

w_col=w%w_pace;
// WriteResult (1, Result_hyp, Result_gamma, w_col, snp_pos, pos_loglr);
WriteResult (1, Result_hyp, Result_gamma, w_col);

gsl_vector *alpha=gsl_vector_alloc (ns_test);
gsl_vector_set_zero (alpha);
//WriteParam (beta_g, alpha, w);
WriteParam (beta_g, alpha, w, snp_pos, pos_loglr);
gsl_vector_free(alpha);

gsl_matrix_free(Result_hyp);
gsl_matrix_free(Result_gamma);

gsl_vector_free(sigma_subvec_old);
gsl_vector_free(sigma_subvec_new);

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
} */


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


/*double BSLMM::ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
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
	
	cHyp_new.h=h;
	cHyp_new.rho=rho;
	return 0.0;
}*/

//specific for two function types
double BSLMM::ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
    double h=cHyp_old.h, rho=cHyp_old.rho_vec[0];
    h_max = 0.99; h_min = 0.01;
	//double d_h=min(0.01, (h_max-h_min)*h_scale);
    double d_h = 0.1;
    
    double rmax = 0.99, rmin = 0.01;
    double rho_h = min(0.01, (rmax - rmin) * rho_scale);
	
	for (size_t i=0; i<repeat; ++i) {
		h=h+(gsl_rng_uniform(gsl_r)-0.5)*d_h;
		if (h<h_min) {h=2*h_min-h;}
		if (h>h_max) {h=2*h_max-h;}
		
		rho=rho+(gsl_rng_uniform(gsl_r)-0.5)*rho_h;
		if (rho<rmin) {rho=2*rmin-rho;}
		if (rho>rmax) {rho=2*rmax-rho;}
	}
    
	cHyp_new.h=h;
	cHyp_new.rho_vec[0]=rho;
    cHyp_new.rho_vec[1]=1.0 - rho;
    
	return 0.0;
}

double BSLMM::ProposeRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
    double rho=cHyp_old.rho_vec[0];
    double rmax = 0.99, rmin = 0.01;
    double rho_h = min(0.01, (rmax - rmin) * rho_scale);
    
    for (size_t i=0; i<repeat; ++i) {
        rho=rho+(gsl_rng_uniform(gsl_r)-0.5)*rho_h;
        if (rho<rmin) {rho=2*rmin-rho;}
        if (rho>rmax) {rho=2*rmax-rho;}
    }
    
    cHyp_new.rho_vec[0]=rho;
    cHyp_new.rho_vec[1]=1.0 - rho;
    
    return 0.0;
}


double BSLMM::CalcPrho(const class HYPBSLMM &cHyp)
{
    //assume dirichlet prior for rho_vec
    double prho = 0.0;
    for (size_t i=0; i < n_type; i++) {
        prho += (double)(mFunc[i] - 1) * log(cHyp.rho_vec[i]);
    }
    return prho;
}

double BSLMM::CalcPh(const class HYPBSLMM &cHyp)
{
    //assume dirichlet prior for rho_vec
    double ph;
    //ph = -log(cHyp.h);
    ph = 4.0 * log(cHyp.h) + 4.0 * log(1.0 - cHyp.h); // Beta(3, 3)
    return ph;
}

double BSLMM::ProposePi (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	double logp_old=cHyp_old.logp, logp_new=cHyp_old.logp;
	double log_ratio = 0.0;
    
	double d_logp=min(0.5, (logp_max-logp_min)*logp_scale);
	
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
    //logp_max = log(0.1);
    vector<double> logtheta_minvec;
    logtheta_minvec.push_back( -log((double)mFunc[0]) );
    logtheta_minvec.push_back( -log((double)mFunc[1]) );
	double d_logp = (logp_max-logp_min)*logp_scale;
	
	for (size_t i=0; i<repeat; ++i) {
        for (size_t j=0; j<n_type; j++) {
            logp_new[j] = logp_old[j] + (gsl_rng_uniform(gsl_r)-0.5)*d_logp;
            if (logp_new[j]<logtheta_minvec[j])
                {logp_new[j]=2*logtheta_minvec[j]-logp_new[j];}
            if (logp_new[j]>logp_max) {logp_new[j]=2*logp_max-logp_new[j];}
            log_ratio+=logp_new[j]-logp_old[j];
            logp_old[j]=logp_new[j];
        }
	}
	cHyp_new.log_theta=logp_new;
	
	return log_ratio;
}

double BSLMM::ProposeSubvar (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	vector<double> subvar_old = cHyp_old.subvar, subvar_new=cHyp_old.subvar;
	double log_ratio = 0.0;
    	
	for (size_t i=0; i<repeat; ++i) {
        for (size_t j=0; j<n_type; j++) {
            subvar_new[j] = gsl_ran_lognormal(gsl_r, log(subvar_old[j]), vscale);
            log_ratio += log(subvar_new[j]) - log(subvar_old[j]);
            subvar_old[j] = subvar_new[j];
        }
	}
	cHyp_new.subvar=subvar_new;
	
	return log_ratio;
}

double BSLMM::ProposeTheta (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, size_t j, gsl_vector *pi_vec_old, gsl_vector *pi_vec_new, const vector<SNPPOS> &snp_pos)
{
    double logp_old = cHyp_old.log_theta[j], logp_new=cHyp_old.log_theta[j];
    double log_ratio = 0.0;
    double d_logp = min(0.5, (logp_max-logp_min)*logp_scale);
    //cout << d_log;
    
    for (size_t i=0; i<repeat; ++i) {
            logp_new = logp_old + (gsl_rng_uniform(gsl_r)-0.5)*d_logp;
            if (logp_new<logp_min) {logp_new=2*logp_min-logp_new;}
            if (logp_new>logp_max) {logp_new=2*logp_max-logp_new;}
            log_ratio+=logp_new-logp_old;
            logp_old=logp_new;
    }
    
    cHyp_new.log_theta[j]=logp_new;
    
    //calculate pi_new_vec
    double dtheta, pi_temp;
    dtheta = exp(cHyp_new.log_theta[j]) - exp(cHyp_old.log_theta[j]);
    
    for (size_t i=0; i < ns_test; i++) {
        if(snp_pos[i].indicator_func[j])
            {
                pi_temp = gsl_vector_get(pi_vec_new, i);
                pi_temp += dtheta * snp_pos[i].weight[j];
                if (pi_temp > 1 || pi_temp < 0) {
                    cerr << "pi_temp = " << pi_temp << endl;
                    exit(-1);
                }
                else {
                    gsl_vector_set(pi_vec_new, i, pi_temp);
                }
            }
    }
    
    return log_ratio;
}

double BSLMM::ProposeSubvar (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat, size_t j, gsl_vector *sigma_vec_new, const vector<size_t> &rank, const vector<SNPPOS> &snp_pos)
{
    double subvar_old = cHyp_old.subvar[j], subvar_new=cHyp_old.subvar[j];
    double log_ratio = 0.0;
    
    for (size_t i=0; i<repeat; ++i) {
            subvar_new = gsl_ran_lognormal(gsl_r, log(subvar_old), vscale);
            log_ratio += log(subvar_new) - log(subvar_old);
            subvar_old = subvar_new;
    }
    cHyp_new.subvar[j]=subvar_new;
    
    size_t order_i;
    double sigma_temp;
   // double d_sigma = cHyp_new.subvar[j] - cHyp_old.subvar[j];
    
    for (size_t i=0; i < rank.size(); i++) {
        order_i = mapRank2Order[rank[i]];
        
        if(snp_pos[order_i].indicator_func[j])
            {
                sigma_temp = CalcSigma(cHyp_new, order_i, snp_pos);
                //sigma_temp = gsl_vector_get(sigma_vec_new, i);
                //sigma_temp += d_sigma * snp_pos[i].weight[j];
          
                if (sigma_temp <= 0) {
                    cerr << "ERROR: sigma_j = " << sigma_temp << "; weight_j = " << snp_pos[order_i].weight[j] << "; old sigma_vec_i = " << gsl_vector_get(sigma_vec_new, i) << endl;
                    exit(-1);
                }
                else gsl_vector_set(sigma_vec_new, i, sigma_temp);
            }
    }
    
    return log_ratio;
}


//JY edit start

void BSLMM::CalcRes(const gsl_matrix *Xgamma, const gsl_vector *z, const gsl_matrix *XtX, const gsl_vector *Xtz, gsl_vector *z_res, const size_t &s_size, const double &ztz){
    
    gsl_matrix_const_view X_gsub=gsl_matrix_const_submatrix(Xgamma, 0, 0, Xgamma->size1, s_size);
    gsl_matrix_const_view XtX_gsub = gsl_matrix_const_submatrix(XtX, 0, 0, s_size, s_size);
    gsl_vector_const_view Xtz_gsub = gsl_vector_const_subvector(Xtz, 0, s_size);
    
    gsl_vector *beta_gamma_hat = gsl_vector_alloc(s_size);
    
    gsl_matrix *XtXtemp = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(XtXtemp, &XtX_gsub.matrix);
    
    double lambda = 0.0;
    for (size_t i=0; i<s_size; ++i) {
        lambda += gsl_matrix_get(XtX, i, i);
    }
    lambda /= (double)s_size;
    lambda *= 0.00001;
   // cout << "labmda = " << lambda << endl;
    gsl_vector_view XtXtemp_diag = gsl_matrix_diagonal(XtXtemp);
    double SSR, R2 ;
    do{
    gsl_vector_add_constant(&XtXtemp_diag.vector, lambda);
    
    if(LapackSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat) != 0)
        EigenSolve(XtXtemp, &Xtz_gsub.vector, beta_gamma_hat);
    //EigenSolve(&XtX_gsub.matrix, &Xtz_gsub.vector, beta_gamma_hat, lambda);
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &X_gsub.matrix, beta_gamma_hat, 0.0, z_res);
    gsl_vector_scale(z_res, -1.0);
    gsl_vector_add(z_res, z);
    
    gsl_blas_ddot(z_res, z_res, &SSR);
    R2 = 1.0 - (SSR / ztz);
    // R2 = R2 - (1 - R2) * s_size / (ni_test - s_size - 1);
   // cout << "R2 = "<< R2 << endl;
   // if(R2 < 0.0) {
      //  cout << "R2 = " << setprecision(6) << R2 << " < 0, "<< endl;
        //cout << "beta_hat estimate from calculating residuals: \n ";
        //PrintVector(beta_gamma_hat);
        //cout << "SSR = " << SSR << "; ztz = " << ztz << "\n";
        //cout << "Set z_res equal to z... " << endl;
      //  gsl_vector_memcpy(z_res, z);
        //NormRes(z_res);
       // }
    }while (R2 < 0.0);
    
    gsl_matrix_free(XtXtemp);
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

//calculate likelihood ratio statistic
double BSLMM::CalcLR(const gsl_vector *z_res, const gsl_vector *x_vec, size_t posj){
    double LR;
    double xtz_res, ztz_res, xtx_vec = XtX_diagvec[posj];
    
    gsl_blas_ddot(z_res, z_res, &ztz_res);
    //gsl_blas_ddot(x_vec, x_vec, &xtx_vec);
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
        rank_j = SNPorder_vec[orderj].second;
        if((orderj >= 0) && (orderj < (long int)ns_test) && (j != win) && (mapRank2in.count(rank_j) == 0)){
            posj = SNPorder_vec[orderj].first;
            getGTgslVec(X, xvec, posj, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
            p_BF[j]=CalcLR(z_res, xvec, posj); //calc loglr
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


bool BSLMM::ColinearTest(uchar ** X, const gsl_matrix * Xtemp, const gsl_matrix * XtX_temp, size_t r_add, size_t s_size)
{
    bool colinear = 0;
    double vreg;
    //double xtx;
    size_t pos = SNPrank_vec[r_add].first;
    double xtx = XtX_diagvec[pos];
    
    gsl_vector *xvec_temp = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec_temp, pos, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
    gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xtemp, 0, 0, Xtemp->size1, s_size);
    gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX_temp, 0, 0, s_size, s_size);
    
    gsl_vector *beta_temp = gsl_vector_alloc(s_size);
    gsl_vector *Xtx_temp = gsl_vector_alloc(s_size);
        
    gsl_blas_dgemv(CblasTrans, 1.0, &Xgamma_sub.matrix, xvec_temp, 0.0, Xtx_temp);
    

    if (LapackSolve(&XtX_sub.matrix, Xtx_temp, beta_temp) !=0) {
        //WriteMatrix(&Xgamma_sub.matrix, "_X");
        //WriteMatrix(&XtX_sub.matrix, "_XtX");
        //WriteVector(Xtx_temp, "_Xtx");
        EigenSolve(&XtX_sub.matrix, Xtx_temp, beta_temp);
        //WriteVector(beta_temp, "_beta");
        
    }
    
    //gsl_blas_ddot(xvec, xvec, &xtx);
    gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
    //cout << "vreg = " << vreg << endl;
    
    double R2 = (vreg / xtx);
    //cout << "R2 = " << R2 << endl;
    if ( R2 >= 0.85 || R2 < 0) {
        colinear = 1;
        /*if (R2 < 0) {
            PrintMatrix(&XtX_sub.matrix, s_size, s_size);
            PrintVector(Xtx_temp);
            PrintVector(beta_temp);
            exit(-1);
        }*/
    }
    
    gsl_vector_free(xvec_temp);
    gsl_vector_free(beta_temp);
    gsl_vector_free(Xtx_temp);
    
    return colinear;
}

/*bool BSLMM::ColinearTest(uchar ** X, gsl_matrix * Xtemp, gsl_matrix * XtX_temp, size_t r_add, size_t s_size)
{
    bool colinear = 0;
    double vreg;
    //double xtx;
    size_t pos = SNPrank_vec[r_add].first;
    double xtx = XtX_diagvec[pos];
    
    gsl_vector *xvec_temp = gsl_vector_alloc(ni_test);
    getGTgslVec(X, xvec_temp, pos, ni_test, ns_test, SNPsd, CompBuffSizeVec, UnCompBufferSize, Compress_Flag);
    gsl_matrix_view Xgamma_sub=gsl_matrix_submatrix (Xtemp, 0, 0, Xtemp->size1, s_size);
    gsl_matrix_view XtX_sub=gsl_matrix_submatrix (XtX_temp, 0, 0, s_size, s_size);

    gsl_vector *beta_temp = gsl_vector_alloc(s_size);
    gsl_vector *Xtx_temp = gsl_vector_alloc(s_size);
    
    gsl_blas_dgemv(CblasTrans, 1.0, &Xgamma_sub.matrix, xvec_temp, 0.0, Xtx_temp);
    
    gsl_matrix *XtXlu = gsl_matrix_alloc(s_size, s_size);
    gsl_matrix_memcpy(XtXlu, &XtX_sub.matrix);
    if (LapackSolve(XtXlu, Xtx_temp, beta_temp) !=0) {
        //WriteMatrix(&Xgamma_sub.matrix, "_X");
        //WriteMatrix(&XtX_sub.matrix, "_XtX");
        //WriteVector(Xtx_temp, "_Xtx");
        
        EigenSolve(&XtX_sub.matrix, Xtx_temp, beta_temp);
        //WriteVector(beta_temp, "_beta");

    }
    gsl_matrix_free(XtXlu);
    
    //gsl_blas_ddot(xvec, xvec, &xtx);
    gsl_blas_ddot(Xtx_temp, beta_temp, &vreg);
    cout << "vreg = " << vreg << endl;
    
    double R2 = (vreg / xtx);
    cout << "R2 = " << R2 << endl;
    if ( R2 >= 0.85 || R2 < 0) {
        colinear = 1;
        if (R2 < 0) {
            PrintMatrix(&XtX_sub.matrix, s_size, s_size);
            PrintVector(Xtx_temp);
            PrintVector(beta_temp);
        }
    }
    
    
    gsl_vector_free(xvec_temp);
    gsl_vector_free(beta_temp);
    gsl_vector_free(Xtx_temp);
    
    return colinear;
}*/


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
       // getGTgslVec(X, &Xgamma_col.vector, pos, Xgamma->size1, ns_test);
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
         
         EigenSolve(Omega, Xty_temp, beta_hat);
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
  //cHyp.pge=1.0;
  return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_vector *Xb, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	
	gsl_blas_ddot (Xb, Xb, &d);
	cHyp.pve=d/(double)ni_test;
	//cHyp.pve/=cHyp.pve+1.0;
	//cHyp.pge=1.0;
	
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

void BSLMM::CreateSnpPosVec(vector<SNPPOS> &snp_pos)
{
    size_t pos;
    string rs;
    string chr;
    long int bp;
    size_t tt=0;
    vector<bool> indicator_func;
    vector<double> weight;
    double weight_i;
    double maf;
    SNPsd.clear();
    
    for (size_t i=0; i < ns_total; ++i){
        if(!indicator_snp[i]) {continue;}
        
        pos=tt;
       // if(tt == pos_loglr[tt].first ) pos = tt;
        //else cout << "error assigning position to snp_pos vector"<< endl;
        
        rs = snpInfo[i].rs_number;
        chr = snpInfo[i].chr;
        bp = snpInfo[i].base_position;
        maf = snpInfo[i].maf;
        indicator_func = snpInfo[i].indicator_func;
        weight = snpInfo[i].weight;
        
        weight_i = snpInfo[i].weight_i;
        // cout << weight_i << ", ";
        SNPPOS snp_temp={pos, rs, chr, bp, maf, indicator_func, weight, weight_i};
        snp_pos.push_back(snp_temp);
        
        SNPsd.push_back(1.0 / sqrt(2.0 * maf * (1.0 - maf)));
        
        tt++;
    }
    snpInfo.clear();
    
}

void BSLMM::InitialMap(const vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos){
    
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

void BSLMM::MHsave(const size_t &t, size_t &w, gsl_matrix *Result_hyp, gsl_matrix *Result_gamma, LModel &model_old, vector<pair<double, double> > &beta_g, const double &mean_z, const vector<pair<size_t, double> > &pos_loglr, const vector<SNPPOS> &snp_pos, double &logPost_old)
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

//if a_mode==13, then run probit model, MCMC with model struct
/*void BSLMM::MCMC (uchar **X_Genotype, gsl_vector *z) {
    
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
    vector<SNPPOS> snp_pos;
    CreateSnpPosVec(snp_pos);
   // printSNPInfo(snp_pos, 10); // print out a few markers' info
	
    InitialMap(pos_loglr, snp_pos);
    double *p_gamma = new double[ns_test];
	CalcPgamma (p_gamma); // calculate discrete distribution for gamma
	CreateGammaProposal(p_gamma);
	
	//initial parameters
    cout << "start initializing mcmc...\n";
    
	//InitialMCMC (X_Genotype, z, model_old, pos_loglr, snp_pos);
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
                cerr << "acceptance percentage = " << accept_percent << " < 0.000001;";
               // exit(1);
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
}*/
