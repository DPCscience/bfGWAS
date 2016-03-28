/*
    Scalable Functional Bayesian Association --- MCMC (SFBA:MCMC)
    Copyright (C) 2016  Jingjing Yang

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

#include "param.h"


typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef short int int16;

using namespace std;

void ProgressBar (string str, double p, double total);

void ProgressBar (string str, double p, double total, double ratio);
std::istream& safeGetline(std::istream& is, std::string& t);

bool ReadFile_snps (const string &file_snps, set<string> &setSnps);

bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo);

bool ReadFile_fam (const string &file_fam, vector<bool> &indicator_pheno, vector<double> &pheno, map<string, int> &PhenoID2Ind, vector<string> & InputSampleID);

bool ReadFile_anno (const string &file_bim, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM);

bool ReadFile_pheno (const string &file_pheno, vector<bool> &indicator_pheno, vector<double > &pheno, vector<string> &InputSampleID);

bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo, size_t &ns_test);

bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, size_t &ns_test);

void ReadFile_kin (const string &file_kin, vector<bool> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G);

bool VCFKin (const string &file_vcf, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin); //** NEED to be rewritten **

bool BimbamKin (const string &file_geno, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);

bool PlinkKin (const string &file_bed, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);

bool ReadFile_geno (const string &file_geno, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar **X, gsl_matrix *K, const bool calc_K, size_t ni_test, size_t ns_test);

bool ReadFile_bed (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar **X, gsl_matrix *K, const bool calc_K, size_t ni_test, size_t ns_test, vector <size_t> &CompBuffSizeVec, bool Compress_Flag);


bool CountFileLines (const string &file_input, size_t &n_lines);

bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, string &GTfield, const map<string, size_t> &PhenoID2Ind, vector<string> &VcfSampleID, vector<size_t> &SampleVcfPos); // first time

bool ReadFile_vcf (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** X, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, string &GTfield, vector<double> &SNPmean, vector <size_t> &CompBuffSizeVec, const vector <size_t> &SampleVcfPos, const map<string, size_t> &PhenoID2Ind, const vector<string> &VcfSampleID, bool Compress_Flag); // second time

bool CreatVcfHash(const string &file_vcf, StringIntHash &sampleID2vcfInd, const string &file_sample);


//for new model
bool ReadFile_anno (const string &file_anno, const string &file_func_code, map<string, int> &mapFunc2Code, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, size_t &n_type, vector<size_t> &mFunc);

void SetMAFCode (const double &maf, string &func_type);

void WriteMatrix(const gsl_matrix * X, const string file_str);
void WriteVector(const gsl_vector * X, const string file_str);

#endif







