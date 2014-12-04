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

#ifndef __IO_H__                
#define __IO_H__


#include <vector>
#include <map>
#include <algorithm>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif

typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef short int int16;

using namespace std;

void ProgressBar (string str, double p, double total);
void ProgressBar (string str, double p, double total, double ratio);
std::istream& safeGetline(std::istream& is, std::string& t);

bool ReadFile_snps (const string &file_snps, set<string> &setSnps);
bool ReadFile_log (const string &file_log, double &pheno_mean);

bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo);
bool ReadFile_fam (const string &file_fam, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, map<string, int> &mapID2num, const vector<size_t> &p_column);

bool ReadFile_cvt (const string &file_cvt, vector<bool> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt);
bool ReadFile_anno (const string &file_bim, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM);
bool ReadFile_pheno (const string &file_pheno, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column);
bool ReadFile_column (const string &file_pheno, vector<bool> &indicator_idv, vector<double> &pheno, const int &p_column);

bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo, size_t &ns_test);
bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test);

void ReadFile_kin (const string &file_kin, vector<bool> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G);
void ReadFile_eigenU (const string &file_u, bool &error, gsl_matrix *U);
void ReadFile_eigenD (const string &file_d, bool &error, gsl_vector *eval); 

bool BimbamKin (const string &file_geno, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);
bool PlinkKin (const string &file_bed, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);

bool ReadFile_geno (const string &file_geno, vector<bool> &indicator_idv, vector<bool> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);
bool ReadFile_bed (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);

bool ReadFile_est (const string &file_est, const vector<size_t> &est_column, map<string, double> &mapRS2est);

bool CountFileLines (const string &file_input, size_t &n_lines);

bool ReadFile_gene (const string &file_gene, vector<double> &vec_read, vector<SNPINFO> &snpInfo, size_t &ng_total);

bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, vector<String> &InputSampleID, StringIntHash &sampleID2vcfInd);

bool ReadFile_vcf (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** UtX, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, vector<String> &InputSampleID, StringIntHash &sampleID2vcfInd);



bool CreatVcfHash(const string &file_vcf, StringIntHash &sampleID2vcfInd);

bool ReadFile_vcf_pheno (const string &file_vcf_pheno, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column, vector<String> &InputSampleID);

bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, const gsl_matrix *W, const gsl_matrix *WtW, const gsl_matrix *WtWi, gsl_vector *Wtx, gsl_vector *WtWiWtx, const vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, const vector<uint> &SampleVcfPos);
bool ReadFile_vcfs (const string &file_vcfs, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** UtX, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, vector<String> &InputSampleID, StringIntHash &sampleID2vcfInd);
#endif







