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
#include <string>
#include <iomanip>
#include <bitset>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"
#include "gzstream.h"
#include "mathfunc.h"
#include "ReadVCF.h"
#include "bslmm.h"
#include "compress.h"

#ifdef FORCE_FLOAT
#include "io_float.h"
#else
#include "io.h"
#endif


using namespace std;



//Print process bar
void ProgressBar (string str, double p, double total)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <50; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%\r"<<flush;
	
	return;
}


//Print process bar (with acceptance ratio)
void ProgressBar (string str, double p, double total, double ratio)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <50; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%    "<<ratio<<"\r"<<flush;
	
	
	return;
}

// in case files are ended with "\r" or "\r\n"
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

//Read snp file
bool ReadFile_snps (const string &file_snps, set<string> &setSnps)
{
	setSnps.clear();

	ifstream infile (file_snps.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open snps file: "<<file_snps<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		setSnps.insert(ch_ptr); 
	}
	
    infile.clear();
	infile.close();
	
	return true;
}


//Read log file
bool ReadFile_log (const string &file_log, double &pheno_mean)
{
	ifstream infile (file_log.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open log file: "<<file_log<<endl; return false;}
	
	string line;
	char *ch_ptr;
	size_t flag=0;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		if (ch_ptr!=NULL && strcmp(ch_ptr, "estimated")==0) {
			ch_ptr=strtok (NULL, " , \t");
			if (ch_ptr!=NULL && strcmp(ch_ptr, "mean")==0) {
				ch_ptr=strtok (NULL, " , \t");
				if (ch_ptr!=NULL && strcmp(ch_ptr, "=")==0) {
					ch_ptr=strtok (NULL, " , \t");
					pheno_mean=atof(ch_ptr);
					flag=1;
				}
			}
		}
		
		if (flag==1) {break;}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}


//Read bimbam annotation file
bool ReadFile_anno (const string &file_anno, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM)
{
	mapRS2chr.clear();
	mapRS2bp.clear();
	
	ifstream infile (file_anno.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening annotation file: "<<file_anno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	long int b_pos;
	string chr;
	double cM;
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		if (strcmp(ch_ptr, "NA")==0) {b_pos=-9;} else {b_pos=atol(ch_ptr);}
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr==NULL || strcmp(ch_ptr, "NA")==0) {chr="-9";} else {chr=ch_ptr;}
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr==NULL || strcmp(ch_ptr, "NA")==0) {cM=-9;} else {cM=atof(ch_ptr);}
		
		mapRS2chr[rs]=chr;
		mapRS2bp[rs]=b_pos;
		mapRS2cM[rs]=cM;
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}

//Read function annotation file
bool ReadFile_anno (const string &file_anno, const string &file_func_code, map<string, int> &mapFunc2Code, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, size_t &n_type)
{
    string line;
    char *pch, *nch;

    //load in unique function codes
    string func_type;
    int func_code;
    
    ifstream infile_code (file_func_code.c_str(), ifstream::in);
    if (!infile_code) {cout<<"error opening annotation file: "<<file_func_code<<endl; return false;}
    while (!safeGetline(infile_code, line).eof()) {
        
        if (line[0] == '#') {
            pch = (char *)line.c_str();
            nch = strchr(pch, '\t');
            n_type = strtol(nch, NULL, 0);
            continue;
        }
        else {
            pch = (char *)line.c_str();
            nch = strchr(pch, '\t');
            func_type.assign(pch, nch-pch);
            func_code = strtol(nch, NULL, 0);
            //cout << func_type << ":" << func_code << endl;
            mapFunc2Code[func_type] = func_code;
        }
    }
    infile_code.close();
    infile_code.clear();
    
    ifstream infile (file_anno.c_str(), ifstream::in);
    if (!infile) {cout<<"error opening annotation file: "<<file_anno<<endl; return false;}
    
    // read function annotation file
    string rs;
    long int b_pos;
    string chr;
    int snp_nfunc;
    size_t snp_i = 0;
    
    while (!safeGetline(infile, line).eof()) {
        if (line[0] == '#') {
            continue;
        }
        else {
          if (!indicator_snp[snp_i]) {
              pch=(char *)line.c_str();
              nch = strchr(pch, '\t');
              rs.assign(pch, nch-pch);
              if (snpInfo[snp_i].rs_number.compare(rs) != 0) {
                  cerr << "annotation file ID dose not match vcf file ID...\n";
                  return false;
              }
              snp_i++;
              continue;
          }
          else{
            pch=(char *)line.c_str();
            nch = strchr(pch, '\t');
            rs.assign(pch, nch-pch);
            if (snpInfo[snp_i].rs_number.compare(rs) != 0) {
                cerr << "annotation file ID dose not match vcf file ID...\n";
                return false;
            }
            
            pch = (nch == NULL) ? NULL : nch+1;
            nch = strchr(pch, '\t');
            chr.assign(pch, nch-pch);
          //  if (pch == NULL || chr.compare("NA") == 0) {
           //     chr = "-9";
           // }
            
            pch = (nch == NULL) ? NULL : nch+1;
            nch = strchr(pch, '\t');
            b_pos = strtol(pch, NULL, 0);
           // if (pch == NULL) {
           //     b_pos = -9;
           // }
              
            pch = (nch == NULL) ? NULL : nch+1;
            snp_nfunc = 0;
            snpInfo[snp_i].indicator_func.assign(n_type, 0);
            //if (snp_i < 5)  cout << rs << ":chr" << chr << ":bp"<< b_pos <<endl;
            while (pch != NULL) {
                nch = strchr(pch, ',');
                if (nch == NULL) func_type.assign(pch);
                else func_type.assign(pch, nch-pch);
                func_code = mapFunc2Code[func_type];
               // cout << func_type << " with code " << func_code << endl;
                if(!snpInfo[snp_i].indicator_func[func_code])
                {
                    snpInfo[snp_i].indicator_func[func_code] = 1;
                    snp_nfunc++;
                }
                pch = (nch == NULL) ? NULL : nch+1;
            }
            
            //if ((snp_nfunc > 0) && (snp_nfunc <= n_type))
              if (snp_nfunc == 1)
              {
                  snpInfo[snp_i].weight_i = 1.0 ;// / (double)snp_nfunc;
                  // CalcWeight(snpInfo[snp_i].indicator_func, snpInfo[snp_i].weight, snpInfo[snp_i].weight_i);
              }
            else if (snp_nfunc == 0) {
                snpInfo[snp_i].weight_i = 0.0;
                indicator_snp[snp_i] = 0;
                cout << "function annotation is NULL \n ";
            }
            else {cerr << "ERROR: snp_nfunc = " <<snp_nfunc<< " ... \n"; exit(-1);}
            snp_i++;
          }
        }
    }
    cout << "n_type = " << n_type << endl;
   // cout << "total snp number = " << snp_i << endl;
    
    infile.close();
    infile.clear();	
    
    return true;
}


//read one column of phenotype
bool ReadFile_column (const string &file_pheno, vector<bool> &indicator_idv, vector<double> &pheno, const int &p_column)
{
	indicator_idv.clear();
	pheno.clear();
	
	igzstream infile (file_pheno.c_str(), igzstream::in);
//	ifstream infile (file_pheno.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_pheno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string id;
	double p;
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		for (int i=0; i<(p_column-1); ++i) {
			ch_ptr=strtok (NULL, " , \t");	
		}		
		if (strcmp(ch_ptr, "NA")==0) {indicator_idv.push_back(0); pheno.push_back(-9);}		//pheno is different from pimass2
		else {p=atof(ch_ptr); indicator_idv.push_back(1); pheno.push_back(p);}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}

//Read VCF phenotype file, p_column=1, 2 ...
bool ReadFile_vcf_pheno (const string &file_vcf_pheno, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column, vector<string> &InputSampleID)
{
	indicator_pheno.clear();
	pheno.clear();
	
   // cout << "open phenotype file ... " << file_vcf_pheno << "\n";
    
	igzstream infile (file_vcf_pheno.c_str(), igzstream::in);
    //	ifstream infile (file_pheno.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_vcf_pheno<<endl; return false;}
    
	string line;
	char *ch_ptr;
    
	string id;
	double p;
	
	vector<double> pheno_row;
	vector<bool> ind_pheno_row;
	
    
	size_t p_max=*max_element(p_column.begin(), p_column.end() );
	map<size_t, size_t> mapP2c;
    
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}
	
    size_t numPheno=0;
	while (!safeGetline(infile, line).eof()) {
        
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
        id=ch_ptr;
		InputSampleID.push_back(id); //load first column as Sample IDs.
        //cout << "id = " << ch_ptr << ":";
        
        ch_ptr=strtok (NULL, " , \t");
        
		size_t i=0;
		while (i<p_max ) {
			if (mapP2c.count(i+1)!=0) {
				if (strcmp(ch_ptr, "NA")==0) {ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;}
                else
                {
                    p=atof(ch_ptr); ind_pheno_row[mapP2c[i+1]]=1;
                    pheno_row[mapP2c[i+1]]=p;
                }
			}
			i++;
			ch_ptr=strtok (NULL, " , \t");
		}
		
		indicator_pheno.push_back(ind_pheno_row);
        //PrintVector(pheno_row);
		pheno.push_back(pheno_row);
        numPheno++;
	}
    cout << "Load numPheno = " << numPheno << "\n";
    
	infile.close();
	infile.clear();
	
	return true;
}





//Read bimbam/VCF phenotype file, p_column=1, 2 ...
bool ReadFile_pheno (const string &file_pheno, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column)
{
	indicator_pheno.clear();
	pheno.clear();
	
    cout << "open phenotype file ... " << file_pheno << "\n";
    
	igzstream infile (file_pheno.c_str(), igzstream::in);
//	ifstream infile (file_pheno.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_pheno<<endl; return false;}

	string line;
	char *ch_ptr;
  
	string id;
	double p;
	
	vector<double> pheno_row;
	vector<bool> ind_pheno_row;
	
    
	size_t p_max=*max_element(p_column.begin(), p_column.end() );
	map<size_t, size_t> mapP2c;
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}	
	
    size_t numPheno=0;
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		
		size_t i=0;
		while (i<p_max ) {			
			if (mapP2c.count(i+1)!=0) {
				if (strcmp(ch_ptr, "NA")==0) {ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;}
                else
                {
                    p=atof(ch_ptr); ind_pheno_row[mapP2c[i+1]]=1;
                    pheno_row[mapP2c[i+1]]=p;
                }
			}
			i++;
			ch_ptr=strtok (NULL, " , \t");	
		}
		
		indicator_pheno.push_back(ind_pheno_row);	
		pheno.push_back(pheno_row);
        numPheno++;
	}
   // cout << "Load numPheno = " << numPheno << "\n";
 
	infile.close();
	infile.clear();	
	
	return true;
}


bool ReadFile_cvt (const string &file_cvt, vector<bool> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt)
{
	indicator_cvt.clear();
	
	ifstream infile (file_cvt.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open covariates file: "<<file_cvt<<endl; return false;}
	
	string line;
	char *ch_ptr;
	double d;	
	
	int flag_na=0;	
	
	while (!safeGetline(infile, line).eof()) {
		vector<double> v_d; flag_na=0;
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		while (ch_ptr!=NULL) {
			if (strcmp(ch_ptr, "NA")==0) {flag_na=1; d=-9;}
			else {d=atof(ch_ptr);}
			
			v_d.push_back(d);
			ch_ptr=strtok (NULL, " , \t");	
		}
		if (flag_na==0) {indicator_cvt.push_back(1);} else {indicator_cvt.push_back(0);} 
		cvt.push_back(v_d);
	}
	
	if (indicator_cvt.empty()) {n_cvt=0;}
	else {
		flag_na=0;
		for (vector<int>::size_type i=0; i<indicator_cvt.size(); ++i) {
			if (indicator_cvt[i]==0) {continue;}
			
			if (flag_na==0) {flag_na=1; n_cvt=cvt[i].size();}
			if (flag_na!=0 && n_cvt!=cvt[i].size()) {cout<<"error! number of covariates in row "<<i<<" do not match other rows."<<endl; return false;}
		}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}



//Read .bim file
bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo)
{
	snpInfo.clear();
	
    cout << "Start reading bim file ...\n";
	ifstream infile (file_bim.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening .bim file: "<<file_bim<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	long int b_pos;
	string chr;
	double cM;
	string major;
	string minor;
    
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		chr=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		cM=atof(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		b_pos=atol(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		major=ch_ptr;
		
        SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, -9, -9, -9, indicator_func_temp, weight_temp, 0.0};
		snpInfo.push_back(sInfo);
	}
	
	infile.close();
	infile.clear();
    cout << "Success reading bim file ...\n";
	return true;
}


//Read .fam file
bool ReadFile_fam (const string &file_fam, vector<vector<bool> > &indicator_pheno, vector<vector<double> > &pheno, map<string, int> &PhenoID2Ind, const vector<size_t> &p_column, vector<string> & InputSampleID)
{
	indicator_pheno.clear();
	pheno.clear();
	PhenoID2Ind.clear();
	
    cout << "Start reading fam file ...\n";
	igzstream infile (file_fam.c_str(), igzstream::in);
	//ifstream infile (file_fam.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening .fam file: "<<file_fam<<endl; return false;}

	string line;
	char *ch_ptr;

	string id;
	size_t c=0;
	double p;

	vector<double> pheno_row;
	vector<bool> ind_pheno_row;
	
	size_t p_max=*max_element(p_column.begin(), p_column.end() );
	map<size_t, size_t> mapP2c;
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}	
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		ch_ptr=strtok (NULL, " \t");
		id=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		
		size_t i=0;
		while (i<p_max ) {
			if (mapP2c.count(i+1)!=0 ) {
				if (strcmp(ch_ptr, "NA")==0) {
					ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;
				} else {
					p=atof(ch_ptr);
					
					if (p==-9) {ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;}
					else {ind_pheno_row[mapP2c[i+1]]=1; pheno_row[mapP2c[i+1]]=p;}
				}
			}
			i++;
			ch_ptr=strtok (NULL, " , \t");	
		}
		
		indicator_pheno.push_back(ind_pheno_row);
		pheno.push_back(pheno_row);				
		InputSampleID.push_back(id);
		PhenoID2Ind[id]=c; c++;
	}
 
	infile.close();
	infile.clear();
    cout << "Success reading fam file ...\n";
	return true;
}

bool CreatVcfHash(const string &file_vcf, StringIntHash &sampleID2vcfInd, const string &file_sample){
        
    VcfFileReader inFile;
    VcfHeader header;
    
    if(!inFile.open(file_vcf.c_str(), header, file_sample.c_str(), NULL, NULL))
    {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(1);
    }
    
    uint numSample = (uint)header.getNumSamples();
    cout << "numSample = " << numSample << endl;
    String sample_name;
	for (size_t i=0; i<numSample; ++i) {
        sample_name = header.getSampleName(i);
        sampleID2vcfInd.Add(sample_name, i);
	}
    cout << "\n create hash sampleID to vcf index success...\n";
    return true;
}

/* void GetVcfPos(const vector<string> &VcfSampleID, const map<string, size_t> &PhenoID2Ind, vector <size_t> &SampleVcfPos)
{
    size_t yidx;
    string sampleid;
    SampleVcfPos.clear();
    
    for (size_t i=0; i < VcfSampleID.size(); i++) {
        sampleid = VcfSampleID[i];
        if (PhenoID2Ind.count(sampleid) == 0) continue;
        else {
            yidx = PhenoID2Ind[sampleid];
            SampleVcfPos.push_back(i);
        }
    }
} */

// Read VCF genotype file, the first time,
bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, string &GTfield, const map<string, size_t> &PhenoID2Ind, vector<string> &VcfSampleID, vector<size_t> &SampleVcfPos)
{
    if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string
    
    indicator_snp.clear();
    snpInfo.clear();
    
    igzstream infile(file_vcf.c_str(), igzstream::in);
    cout << "open vcf file ...\n";
    if(!infile) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(-1);
    }
    
    gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
    gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
    gsl_vector *Wtx=gsl_vector_alloc (W->size2);
    gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
    gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
    int sig;
    LUDecomp (WtW, pmt, &sig);
    LUInvert (WtW, pmt, WtWi);
    double v_x, v_w;
    
    string rs; long int b_pos = 0; string chr;
    string major; string minor; double cM=-9;
    string s;
    
    double maf, geno, geno_old;
    size_t n_miss, n_0, n_1, n_2;
    int flag_poly;
    
    size_t ni_total = indicator_idv.size();
    size_t c_idv=0; //count individual number in each record
    //cout << "ni_total = indicator_idv.size() =" << indicator_idv.size() << endl;
    
    ns_test=0; // variable defined in param.h
    gsl_vector *genotype = gsl_vector_alloc(W->size1);
    // cout << "genotype size = " << W->size1 << "\n" ;
    vector<bool> genotype_miss(ni_test, 0);
    
    char *pch, *p, *nch=NULL, *n;
    size_t tab_count;
    int GTpos=0, k=0;
    
    string line;
    VcfSampleID.clear();
    SampleVcfPos.clear(); // with length = ni_total
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
    
    // cout << "start reading record ... \n";
  while(!safeGetline(infile, line).eof()) {
        if (line[0] == '#') {
           if (strncmp(line.c_str(), "#CHROM", 6) == 0) {
               pch= (char *)line.c_str();
             //parse for individual IDs, save VCFsampleID, create SampleVcfPos
               for (tab_count=0; pch != NULL; tab_count++) {
                   nch=strchr(pch, '\t'); //point to the position of next '\t'
                   if (tab_count>8) {
                       if (nch == NULL) { s.assign( pch );}
                       else s.assign( pch, nch-pch );
                       VcfSampleID.push_back(s);
                       if (PhenoID2Ind.count(s)>0) {
                           SampleVcfPos.push_back(tab_count); //record tab_position
                       }
                   }
                   pch = (nch == NULL) ? NULL : nch+1;
               }
               //cout << "Parse for VCF sample IDs \n";
            }
            continue;
        }
    else{
        c_idv=0; n_0=0; n_1=0; n_2=0;
        maf=0; n_miss=0; flag_poly=0; geno_old=-9;
        
        pch= (char *)line.c_str();
        
        for (tab_count=0; pch != NULL; tab_count++) {
            nch=strchr(pch, '\t'); //point to the position of next '\t'
            if (tab_count<5) {
                if (nch == NULL) { s.assign( pch );}
                else s.assign( pch, nch-pch ); // field string s
                
                switch (tab_count) {
                    case 0:
                        chr = s; break;
                    case 1:
                        b_pos=atol(s.c_str()); break;
                    case 2:
                        rs = s; break;
                    case 3:
                        major = s; break;
                    case 4:
                        minor = s; break;
                    default:
                        break;
                }
                if (setSnps.size()!=0 && setSnps.count(rs)==0) {
                    indicator_snp.push_back(0);
                    continue;
                }
            }
            
            else if ((tab_count == 6) && (pch[0] == 'F')){
                SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
                snpInfo.push_back(sInfo); //save marker information
                indicator_snp.push_back(0);
                continue;
                //failed filter, continue with next record
            }
            
            else if ((tab_count == 8) && (c_idv == 0))
            {
                // parse FORMAT field
                if (pch[0] == GTfield[0] && pch[1] == GTfield[1] && ((nch==pch+2)||pch[2]==':') ) {
                    GTpos=0; //GT start in the first position
                }
                else if (nch == NULL){ cerr << "VCF has FORMAT field but dose not have any genotype\n";}
                else{
                    k=0; //index of key characters
                    GTpos=0;
                    p=pch;
                    while (p<nch) {
                        if (*p == ':') {
                            if (k >= lkey) {
                                break;
                            }
                            else {
                                ++GTpos;
                                k=0;
                            }
                        }
                        else {
                            if (GTfield[k] == *p) {
                                ++k;
                            }
                            else { k=0; }
                        }
                     ++p;
                    }
                    if ((p==nch) && (k != lkey)) {
                        cerr << "Cannot find" << GTfield << "at marker" << chr << ":" << b_pos << endl;
                        exit(-1);
                    }
                }
            }
            else if ( tab_count == SampleVcfPos[c_idv] )
                {
                  if ( !indicator_idv[c_idv] ) {
                      c_idv++; continue;
                  }
                  else{
                    p = pch; // make p reach to the key index
                    if (GTpos>0) {
                        for (int i=0; (i<GTpos) && (p!=NULL); ++i) {
                            n = strchr(p, ':');
                            p = (n == NULL) ? NULL : n+1;
                        }
                    }
                    
                    if ((p==NULL) || ((p[0]=='.') && (p[2]=='.')) ) {
                        geno = -9;//missing
                        genotype_miss[c_idv]=1; n_miss++; c_idv++; continue;
                    }
                    else if ( (p[1] == '/') || (p[1] == '|') ) {
                        //read bi-allelic GT
                        if (p[0]=='.') {
                            geno = (double)(p[2] -'0');
                        }
                        else if (p[2]=='.') {
                            geno = (double)(p[0] -'0');
                        }
                        else geno = (double)((p[0] - '0') + (p[2]- '0'));
                    }
                    else {
                        //read dosage data
                        geno = strtod(p, NULL);
                    }
                    if (geno>=0 && geno<=0.5) {n_0++;}
                    if (geno>0.5 && geno<1.5) {n_1++;}
                    if (geno>=1.5 && geno<=2.0) {n_2++;}
                    gsl_vector_set (genotype, c_idv, geno);
                    if (flag_poly==0) {geno_old=geno; flag_poly=2;}
                    if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
                    maf+=geno;
                    c_idv++;
                  }
                }
            pch = (nch == NULL) ? NULL : nch+1;
        }
        if (c_idv != ni_total) {
            cerr << "record sample number " << c_idv << " dose not equal to ni_total " << ni_total << "\n";
            exit(-1);
        }
        SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
        snpInfo.push_back(sInfo); //save marker information
        
        maf/=2.0*(double)(ni_test-n_miss);
        //cout << "maf = " << maf << "\n";
        
        if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
        //cout << "pass missness criteron...\n";
        
        if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}
        
        if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
        //cout << "pass maf criteron...\n";
        
        if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
        // cout << "pass poly criteron...\n";
        
        if (hwe_level!=0) {
            if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
        }
        //  cout << "pass hwe criteron...\n";
        
        //filter SNP if it is correlated with W
        for (size_t i=0; i<ni_test; ++i) {
            if (genotype_miss[i]) {
                geno=maf*2.0;
                gsl_vector_set (genotype, i, geno);
            }
        }
        
        gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
        gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
        gsl_blas_ddot (genotype, genotype, &v_x);
        gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
        
        if (v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}
        
        indicator_snp.push_back(1);
        ns_test++;
        //if (ns_test < 5) {sInfo.printMarker(); PrintVector(genotype, 10);}
        }
    }
    //cout << "genotype vector:\n";
    // PrintVector(genotype, 10);
    cout << "VCF tab_count = " << tab_count << endl;
    cout << "ns_test = " << ns_test ;
    cout << "vcf read first time success ... \n";
    // cout << "ns_test = " << ns_test << "indicator_snp.size = " << indicator_snp.size()<<"\n";
    
    gsl_vector_free (genotype);
    gsl_matrix_free (WtW);
    gsl_matrix_free (WtWi);
    gsl_vector_free (Wtx);
    gsl_vector_free (WtWiWtx);
    gsl_permutation_free (pmt);
    
    infile.clear();
    infile.close();
    //inFile.clear();
    return true;
}


/* bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, vector<String> &InputSampleID, StringIntHash &sampleID2vcfInd, const string &file_sample)
{
	
    indicator_snp.clear();
    snpInfo.clear();
    
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    
    // Set to only store the GT genotype field.
    VcfRecordGenotype::addStoreField("GT");
    VcfRecordGenotype::addStoreField("EC");
    
    cout << "open vcf file ...\n";
    if(!inFile.open(file_vcf.c_str(), header, file_sample.c_str(), NULL, NULL)) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(1);
    }
    
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);
	
	double v_x, v_w;
	size_t c_idv=0;
    
	string rs;
	long int b_pos;
	string chr;
	string major;
	string minor;
    double cM=-9;
    
	double maf, geno, geno_old;
	size_t n_miss;
	size_t n_0, n_1, n_2;
	int16 flag_poly;
    
    size_t ni_total = indicator_idv.size();
    //cout << "ni_total = indicator_idv.size() =" << indicator_idv.size() << endl;
    vector<uint> SampleVcfPos;
    uint vcfpos;
    
    cout << "SampleVcfPos: " ;
    for (int i=0; i<(int)ni_total; i++) {
        vcfpos = (uint)sampleID2vcfInd.Integer(InputSampleID[i]);
        SampleVcfPos.push_back(vcfpos);
    }
    
    ns_test=0;
    genMarker temp_genMarker;
    
    gsl_vector *genotype = gsl_vector_alloc(W->size1);
   // cout << "genotype size = " << W->size1 << "\n" ;
    vector<bool> genotype_miss(ni_test, 0);

   // cout << "start reading record ... \n";
	while(inFile.readRecord(record)) {
        
        temp_genMarker.iniRecord(record);
        chr = temp_genMarker.chr;
        rs = temp_genMarker.rs;
        b_pos = temp_genMarker.bp;
        minor = temp_genMarker.Alt;
        major = temp_genMarker.Ref;
		
		if (setSnps.size()!=0 && setSnps.count(rs)==0) {
			indicator_snp.push_back(0);
			continue;
		}
        
		maf=0; n_miss=0; flag_poly=0; geno_old=-9;
		n_0=0; n_1=0; n_2=0;
		c_idv=0;
        
        for (size_t i=0; i<ni_total; ++i) {
            
            if (!indicator_idv[i]) {continue;}
            
            geno = getDoubleDosageFromRecord(record, SampleVcfPos[i]);
            if (geno == -9) {genotype_miss[c_idv]=1; n_miss++; c_idv++; continue;}
            
            if (geno>=0 && geno<=0.5) {n_0++;}
            if (geno>0.5 && geno<1.5) {n_1++;}
            if (geno>=1.5 && geno<=2.0) {n_2++;}
            gsl_vector_set (genotype, c_idv, geno);
            
            if (flag_poly==0) {geno_old=geno; flag_poly=2;}
            if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
            
            maf+=geno;
            c_idv++;
            
        }
        
        //PrintVector(genotype);
        
		maf/=2.0*(double)(ni_test-n_miss);
        //cout << "maf = " << maf << "\n";
		
		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf};
		snpInfo.push_back(sInfo);
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
        //cout << "pass missness criteron...\n";
        
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}
        
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
       //cout << "pass maf criteron...\n";
		
		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
       // cout << "pass poly criteron...\n";
		
		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
      //  cout << "pass hwe criteron...\n";
		
		//filter SNP if it is correlated with W
		for (size_t i=0; i<ni_test; ++i) {
			if (genotype_miss[i]) {
                geno=maf*2.0;
                gsl_vector_set (genotype, i, geno);
            }
		}
		
		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
		
		if (v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}
		
		indicator_snp.push_back(1);
		ns_test++;
       // cout << "ns_test = " << ns_test ;
	}
    //cout << "genotype vector:\n";
   // PrintVector(genotype, 10);
    cout << "vcf read first time success ... \n";
   // cout << "ns_test = " << ns_test << "indicator_snp.size = " << indicator_snp.size()<<"\n";
	
	gsl_vector_free (genotype);
	gsl_matrix_free (WtW);
	gsl_matrix_free (WtWi);
	gsl_vector_free (Wtx);
	gsl_vector_free (WtWiWtx);
	gsl_permutation_free (pmt);
    
    
    inFile.close();
    //inFile.clear();
	return true;
}
*/

//read VCF file given multiple file names
// Read VCF genotype file, the first time,
bool ReadFile_vcf (const string &file_vcf, const set<string> &setSnps, const gsl_matrix *W, const gsl_matrix *WtW, const gsl_matrix *WtWi, gsl_vector *Wtx, gsl_vector *WtWiWtx, const vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, vector<SNPINFO> &snpInfo, size_t &ns_test, size_t &ni_test, const vector<uint> &SampleVcfPos, const string &file_sample)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype::addStoreField("GT");
    VcfRecordGenotype::addStoreField("EC");
    if(!inFile.open(file_vcf.c_str(), header, file_sample.c_str(), NULL, NULL)) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(1);
    }
    
    double v_x, v_w;
    size_t c_idv=0;
    
    string rs;
    long int b_pos;
    string chr;
    string major;
    string minor;
    double cM=-9;
    
    double maf, geno, geno_old;
    size_t n_miss;
    size_t n_0, n_1, n_2;
    int16 flag_poly;
    
    genMarker temp_genMarker;
    
    gsl_vector *genotype = gsl_vector_alloc(W->size1);
    vector<bool> genotype_miss(ni_test, 0);
    
    // cout << "start reading record ... \n";
    while(inFile.readRecord(record)) {
        
        temp_genMarker.iniRecord(record);
        chr = temp_genMarker.chr;
        rs = temp_genMarker.rs;
        b_pos = temp_genMarker.bp;
        minor = temp_genMarker.Alt;
        major = temp_genMarker.Ref;
        
        if (setSnps.size()!=0 && setSnps.count(rs)==0) {
            indicator_snp.push_back(0);
            continue;
        }
        
        maf=0; n_miss=0; flag_poly=0; geno_old=-9;
        n_0=0; n_1=0; n_2=0;
        c_idv=0;
        vector<bool> indicator_func_temp;
        vector<double> weight_temp;
        
        for (size_t i=0; i<indicator_idv.size(); ++i) {
            
            if (!indicator_idv[i]) {continue;}
            
            geno = getDoubleDosageFromRecord(record, SampleVcfPos[i]);
            if (geno == -9.0) {genotype_miss[c_idv]=1; n_miss++; c_idv++; continue;}
            
            if (geno>=0 && geno<=0.5) {n_0++;}
            if (geno>0.5 && geno<1.5) {n_1++;}
            if (geno>=1.5 && geno<=2.0) {n_2++;}
            gsl_vector_set (genotype, c_idv, geno);
            
            if (flag_poly==0) {geno_old=geno; flag_poly=2;}
            if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
            
            maf+=geno;
            c_idv++;
        }
        
        maf/=2.0*(double)(ni_test-n_miss);
        
        SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
        snpInfo.push_back(sInfo);
        
        if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
        
        if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
        
        if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
        
        if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

        if (hwe_level!=0) {
            if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
        }
        
        //filter SNP if it is correlated with W
        for (size_t i=0; i<ni_test; ++i) {
            if (genotype_miss[i]) {
                geno=maf*2.0;
                gsl_vector_set (genotype, i, geno);
            }
        }
        
        gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
        gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
        gsl_blas_ddot (genotype, genotype, &v_x);
        gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
        if (v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}
        
        indicator_snp.push_back(1);
        ns_test++;
    }
    cout << "ns_test = " << ns_test << "; indicator_snp.size = " << indicator_snp.size()<<"\n";
    gsl_vector_free (genotype);
    inFile.close();
    return true;
}



//Read bimbam mean genotype file, the first time, to obtain #SNPs for analysis (ns_test) and total #SNP (ns_total)
bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo, size_t &ns_test) {
    
	indicator_snp.clear();
	snpInfo.clear();
	
	igzstream infile (file_geno.c_str(), igzstream::in);
//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}

	gsl_vector *genotype=gsl_vector_alloc (W->size1);
	gsl_vector *genotype_miss=gsl_vector_alloc (W->size1);
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;	
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);
	
	double v_x, v_w;
	size_t c_idv=0;
	
	string line;
	char *ch_ptr;
		
	string rs;
	long int b_pos;
	string chr;
    //string chr_interest = "1 5 8 10 15 16 18";
	string major;
	string minor;
	double cM;
  
	double maf, geno, geno_old;
	size_t n_miss;
	size_t n_0, n_1, n_2;
	int flag_poly;
	
	size_t ni_total=indicator_idv.size();
	size_t ni_test=0;
	for (size_t i=0; i<ni_total; ++i) {
		ni_test+=indicator_idv[i];
	}
	ns_test=0;
	
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
    
	while (!safeGetline(infile, line).eof()) {		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		major=ch_ptr;
		
		if (setSnps.size()!=0 && setSnps.count(rs)==0) {
			SNPINFO sInfo={"-9", rs, -9, -9, minor, major, -9, -9, -9, indicator_func_temp, weight_temp, 0.0};
			snpInfo.push_back(sInfo);
			indicator_snp.push_back(0);
			continue;
		}
				
		if (mapRS2bp.count(rs)==0) {chr="-9"; b_pos=-9;cM=-9;}
		else {b_pos=mapRS2bp[rs]; chr=mapRS2chr[rs]; cM=mapRS2cM[rs];}		
				
		maf=0; n_miss=0; flag_poly=0; geno_old=-9;
		n_0=0; n_1=0; n_2=0; c_idv=0;
        
        gsl_vector_set_zero (genotype_miss);
		for (size_t i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}		

			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++; c_idv++; continue;}
			
			geno=atof(ch_ptr);
			if (geno>=0 && geno<=0.5) {n_0++;}
			if (geno>0.5 && geno<1.5) {n_1++;}
			if (geno>=1.5 && geno<=2.0) {n_2++;}
			
			gsl_vector_set (genotype, c_idv, geno); 
			
//			if (geno<0) {n_miss++; continue;}
			
			if (flag_poly==0) {geno_old=geno; flag_poly=2;}
			if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
			
			maf+=geno;
			
			c_idv++;
		}
		maf/=2.0*(double)(ni_test-n_miss);	
		
		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, (int)n_miss, (double)n_miss/(double)ni_test, maf, indicator_func_temp, weight_temp, 0.0};
		snpInfo.push_back(sInfo);
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
		
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
		
		//filter SNP if it is correlated with W
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get(genotype_miss, i)==1)
            {
                geno=maf*2.0; gsl_vector_set (genotype, i, geno);
            }
		}
		
		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
		
		if (v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}
		
		indicator_snp.push_back(1); 
		ns_test++;
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	gsl_matrix_free (WtW);
	gsl_matrix_free (WtWi);
	gsl_vector_free (Wtx);
	gsl_vector_free (WtWiWtx);
	gsl_permutation_free (pmt);
	
	infile.close();
	infile.clear();	
	
	return true;
}




//Read bed file, the first time
bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test)
{
	indicator_snp.clear();
	size_t ns_total=snpInfo.size();
	
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
    
	gsl_vector *genotype=gsl_vector_alloc (W->size1);
	gsl_vector *genotype_miss=gsl_vector_alloc (W->size1);
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);
	
	double v_x, v_w, geno, geno_old;
	size_t c_idv=0;
	
	char ch[1];
	bitset<8> b;
  	
	size_t ni_total=indicator_idv.size();
	size_t ni_test=0;
	for (size_t i=0; i<ni_total; ++i) {
		ni_test+=indicator_idv[i];
	}
	ns_test=0;
	
	//calculate n_bit and c, the number of bit for each snp
	size_t n_bit;
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}
    
	//ignore the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	double maf, vtx;
	size_t n_miss;
	size_t n_0, n_1, n_2, c;
    int flag_poly;
	
	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		if (setSnps.size()!=0 && setSnps.count(snpInfo[t].rs_number)==0) {
			snpInfo[t].n_miss=-9;
			snpInfo[t].missingness=-9;
			snpInfo[t].maf=-9;
			indicator_snp.push_back(0);
			continue;
		}
        
		//read genotypes
		c=0; maf=0.0; n_miss=0; n_0=0; n_1=0; n_2=0;
        flag_poly=0; geno_old=-9;
		c_idv=0; gsl_vector_set_zero (genotype_miss);
		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c==ni_total) {break;}
				if (indicator_idv[c]==0) {c++; continue;}
				c++;
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); maf+=2.0; n_2++;}
					else {gsl_vector_set(genotype, c_idv, 1.0); maf+=1.0; n_1++;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); maf+=0.0; n_0++;}
					else {gsl_vector_set(genotype_miss, c_idv, 1.0); n_miss++; }
				}
                geno = gsl_vector_get(genotype, c_idv);
                if (flag_poly==0) {geno_old=geno; flag_poly=2;}
                if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
                
				c_idv++;
			}
		}
		maf/=2.0*(double)(ni_test-n_miss);
		
		snpInfo[t].n_miss=n_miss;
		snpInfo[t].missingness=(double)n_miss/(double)ni_test;
		snpInfo[t].maf=maf;
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

        if (flag_poly!=1) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
        
		
		//filter SNP if it is correlated with W
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}
		}
		
        //JY add
        gsl_blas_ddot(genotype, genotype, &vtx);
        if(vtx < 0.00000001)
        {cout << "snp has x'x = " << setprecision(9) << vtx;}
		//JY add
        
		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
		if (v_w/v_x > r2_level) {indicator_snp.push_back(0); continue;}
		
		indicator_snp.push_back(1);
		ns_test++;
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	gsl_matrix_free (WtW);
	gsl_matrix_free (WtWi);
	gsl_vector_free (Wtx);
	gsl_vector_free (WtWiWtx);
	gsl_permutation_free (pmt);
    
	infile.close();
	infile.clear();
	
	return true;
}

      
//Read bed file, the first time
/* bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, const gsl_matrix *W, vector<bool> &indicator_idv, vector<bool> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test)
{
	indicator_snp.clear();
	size_t ns_total=snpInfo.size();
	
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}

	gsl_vector *genotype=gsl_vector_alloc (W->size1);
	gsl_vector *genotype_miss=gsl_vector_all
*/


void ReadFile_kin (const string &file_kin, vector<bool> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G)
{
	igzstream infile (file_kin.c_str(), igzstream::in);
    //	ifstream infile (file_kin.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open kinship file: "<<file_kin<<endl; error=true; return;}
	
	size_t ni_total=indicator_idv.size();
	
	gsl_matrix_set_zero (G);
	
	string line;
	char *ch_ptr;
	double d;
	
	if (k_mode==1) {
		size_t i_test=0, i_total=0, j_test=0, j_total=0;
		while (getline(infile, line)) {
			if (i_total==ni_total) {cout<<"error! number of rows in the kinship file is larger than the number of phentypes."<<endl; error=true;}
			
			if (indicator_idv[i_total]==0) {i_total++; continue;}
			
			j_total=0; j_test=0;
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			while (ch_ptr!=NULL) {
				if (j_total==ni_total) {cout<<"error! number of columns in the kinship file is larger than the number of phentypes for row = "<<i_total<<endl; error=true;}
				
				d=atof(ch_ptr);
				if (indicator_idv[j_total]==1) {gsl_matrix_set (G, i_test, j_test, d); j_test++;}
				j_total++;
				
				ch_ptr=strtok (NULL, " , \t");
			}
			if (j_total!=ni_total) {cout<<"error! number of columns in the kinship file do not match the number of phentypes for row = "<<i_total<<endl; error=true;}
			i_total++; i_test++;
		}
		if (i_total!=ni_total) {cout<<"error! number of rows in the kinship file do not match the number of phentypes."<<endl; error=true;}
	}
	else {
		map<size_t, size_t> mapID2ID;
		size_t c=0;
		for (size_t i=0; i<indicator_idv.size(); i++) {
			if (indicator_idv[i]==1) {mapID2ID[i]=c; c++;}
		}
		
		string id1, id2;
		double Cov_d;
		size_t n_id1, n_id2;
		
		while (getline(infile, line)) {
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			id1=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			id2=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			d=atof(ch_ptr);
			if (mapID2num.count(id1)==0 || mapID2num.count(id2)==0) {continue;}
			if (indicator_idv[mapID2num[id1]]==0 || indicator_idv[mapID2num[id2]]==0) {continue;}
			
			n_id1=mapID2ID[mapID2num[id1]];
			n_id2=mapID2ID[mapID2num[id2]];
			
			Cov_d=gsl_matrix_get(G, n_id1, n_id2);
			if (Cov_d!=0 && Cov_d!=d) {cout<<"error! redundant and unequal terms in the kinship file, for id1 = "<<id1<<" and id2 = "<<id2<<endl;}
			else {
				gsl_matrix_set(G, n_id1, n_id2, d);
				gsl_matrix_set(G, n_id2, n_id1, d);
			}
		}
	}
	
	infile.close();
	infile.clear();
	
	return;
}






void ReadFile_eigenU (const string &file_ku, bool &error, gsl_matrix *U)
{
	igzstream infile (file_ku.c_str(), igzstream::in);
    //	ifstream infile (file_ku.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open the U file: "<<file_ku<<endl; error=true; return;}
	
	size_t n_row=U->size1, n_col=U->size2, i_row=0, i_col=0;
	
	gsl_matrix_set_zero (U);
	
	string line;
	char *ch_ptr;
	double d;
	
	while (getline(infile, line)) {
		if (i_row==n_row) {cout<<"error! number of rows in the U file is larger than expected."<<endl; error=true;}
        
		i_col=0;
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		while (ch_ptr!=NULL) {
			if (i_col==n_col) {cout<<"error! number of columns in the U file is larger than expected, for row = "<<i_row<<endl; error=true;}
			
			d=atof(ch_ptr);
			gsl_matrix_set (U, i_row, i_col, d);
			i_col++;
			
			ch_ptr=strtok (NULL, " , \t");
		}
		
		i_row++;
	}
    
	infile.close();
	infile.clear();
	
	return;
}




void ReadFile_eigenD (const string &file_kd, bool &error, gsl_vector *eval)
{
	igzstream infile (file_kd.c_str(), igzstream::in);
    //	ifstream infile (file_kd.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open the D file: "<<file_kd<<endl; error=true; return;}
	
	size_t n_row=eval->size, i_row=0;
	
	gsl_vector_set_zero (eval);
	
	string line;
	char *ch_ptr;
	double d;
	
	while (getline(infile, line)) {
		if (i_row==n_row) {cout<<"error! number of rows in the D file is larger than expected."<<endl; error=true;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		d=atof(ch_ptr);
		
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr!=NULL) {cout<<"error! number of columns in the D file is larger than expected, for row = "<<i_row<<endl; error=true;}
		
		gsl_vector_set (eval, i_row, d);
		
		i_row++;
	}
	
	infile.close();
	infile.clear();
	
	return;
}



//read bimbam mean genotype file and calculate kinship matrix
bool BimbamKin (const string &file_geno, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	//ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	size_t n_miss;
	double d, geno_mean, geno_var;
	
	size_t ni_total=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_total);
	gsl_vector *geno_miss=gsl_vector_alloc (ni_total);
    
	size_t ns_test=0;
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		!safeGetline(infile, line).eof();
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		geno_mean=0.0; n_miss=0; geno_var=0.0;
		gsl_vector_set_all(geno_miss, 0);
		for (size_t i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set(geno_miss, i, 0); n_miss++;}
			else {
				d=atof(ch_ptr);
				gsl_vector_set (geno, i, d);
				gsl_vector_set (geno_miss, i, 1);
				geno_mean+=d;
				geno_var+=d*d;
			}
		}
		
		geno_mean/=(double)(ni_total-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_total;
		geno_var-=geno_mean*geno_mean;
        //		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_total; ++i) {
			if (gsl_vector_get (geno_miss, i)==0) {gsl_vector_set(geno, i, geno_mean);}
		}
		
		gsl_vector_add_constant (geno, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);}
			else {cout<<"Unknown kinship mode."<<endl;}
		}
		
		ns_test++;
    }
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_total; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno);
	gsl_vector_free (geno_miss);
	
	infile.close();
	infile.clear();
	
	return true;
}







bool PlinkKin (const string &file_bed, vector<bool> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
    
	char ch[1];
	bitset<8> b;
	
	size_t n_miss, ci_total;
	double d, geno_mean, geno_var;
	
	size_t ni_total=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_total);
    
	size_t ns_test=0;
	int n_bit;
	
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }
    
	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		geno_mean=0.0;	n_miss=0; ci_total=0; geno_var=0.0;
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==ni_total) {break;}
                
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(geno, ci_total, 2.0); geno_mean+=2.0; geno_var+=4.0; }
					else {gsl_vector_set(geno, ci_total, 1.0); geno_mean+=1.0; geno_var+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(geno, ci_total, 0.0); }
					else {gsl_vector_set(geno, ci_total, -9.0); n_miss++; }
				}
                
				ci_total++;
			}
		}
        
		geno_mean/=(double)(ni_total-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_total;
		geno_var-=geno_mean*geno_mean;
        //		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_total; ++i) {
			d=gsl_vector_get(geno,i);
			if (d==-9.0) {gsl_vector_set(geno, i, geno_mean);}
		}
		
		gsl_vector_add_constant (geno, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);}
			else {cout<<"Unknown kinship mode."<<endl;}
		}
		
		ns_test++;
    }
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_total; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno);
	
	infile.close();
	infile.clear();
	
	return true;
}


//Read VCF genotype file, the second time, recode genotype and calculate K
bool ReadFile_vcf (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** UtX, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, string &GTfield, vector <size_t> &CompBuffSizeVec, const vector <size_t> &SampleVcfPos)
{
    if (GTfield.empty()) {
        GTfield = "GT"; //defalt load GT Data
    }
    int lkey = GTfield.size(); //length of the field-key string
    
    //size_t ni_total = indicator_idv.size();
    //size_t ns_total = indicator_snp.size();
    
    // Open the VCF file.
    igzstream infile(file_vcf.c_str(), igzstream::in);
    cout << "open vcf file second time ...\n";
    if(!infile) {
        std::cerr << "Unable to open " << file_vcf << "\n";
        exit(-1);
    }
    
    if (calc_K==true) {gsl_matrix_set_zero (K);}
    
    gsl_vector *genotype=gsl_vector_alloc (ni_test);
    uchar *geno_uchar = new uchar[ni_test];
    size_t sourceBufferSize = (ni_test) * sizeof(uchar);
    
    const size_t BufferSize = (size_t)(compressBound(sourceBufferSize));
    uchar * TempCompBuffer = (uchar*)malloc(BufferSize);
    uchar * TempBuffer = (uchar*)malloc(sourceBufferSize);
    
    size_t compressedBufferSize = BufferSize;
    //cout << "Source Buffer Size = " << sourceBufferSize << "; Comp Buffer Bound = " << BufferSize  << endl;
    
    CompBuffSizeVec.clear();
    
    double geno, geno_mean, vtx;
    size_t n_miss, c_idv=0, c_snp=0, ctest_snp = 0, ctest_idv=0;
    int result;
    
    char *pch, *p, *nch=NULL, *n;
    size_t tab_count;
    int GTpos=0, k=0;
    string line;
    
    while(!safeGetline(infile, line).eof())
    {
        if (line[0] == '#') {
            continue; //skip header
        }
        else {
            if (!indicator_snp[c_snp]) {c_snp++; continue;}
            c_idv=0; //increase to the total individuals ni_total
            ctest_idv=0; // increase to the total analyzed individuals
            geno_mean=0.0; n_miss=0;
            vector<bool> genotype_miss(ni_test, 0);
        
            pch= (char *)line.c_str();
            for (tab_count=0; pch != NULL; tab_count++) {
                nch=strchr(pch, '\t'); //point to the position of next '\t'
                if ((tab_count == 8) && (c_idv == 0))
                {
                    // parse FORMAT field
                    if (pch[0] == GTfield[0] && pch[1] == GTfield[1] && ((nch==pch+2)||pch[2]==':') ) {
                        GTpos=0; //GT start in the first position
                    }
                    else if (nch == NULL){ cerr << "VCF has FORMAT field but dose not have any genotype\n";}
                    else{
                        k=0; //index of key characters
                        GTpos=0;
                        p=pch;
                        while (p<nch) {
                            if (*p == ':') {
                                if (k >= lkey) {
                                    break;
                                }
                                else {
                                    ++GTpos;
                                    k=0;
                                }
                            }
                            else {
                                if (GTfield[k] == *p) {
                                    ++k;
                                }
                                else { k=0; }
                            }
                            ++p;
                        }
                        if ((p==nch) && (k != lkey)) {
                            cerr << "Cannot find" << GTfield << endl;
                            exit(-1);
                        }
                    }
                }
                else if ( tab_count == SampleVcfPos[c_idv] )
                {
                    if ( !indicator_idv[c_idv] ) {
                        c_idv++; continue;
                    }
                    else{
                        p = pch; // make p reach to the key index
                        if (GTpos>0) {
                            for (int i=0; (i<GTpos) && (p!=NULL); ++i) {
                                n = strchr(p, ':');
                                p = (n == NULL) ? NULL : n+1;
                            }
                        }
                        
                        if ((p==NULL) || ((p[0]=='.') && (p[2]=='.')) ) {
                            geno = -9;//missing
                            genotype_miss[ctest_idv]=1; n_miss++; c_idv++; continue;
                        }
                        else if ( (p[1] == '/') || (p[1] == '|') ) {
                            //read bi-allelic GT
                            if (p[0]=='.') {
                                geno = (double)(p[2] -'0');
                            }
                            else if (p[2]=='.') {
                                geno = (double)(p[0] -'0');
                            }
                            else geno = (double)((p[0] - '0') + (p[2]- '0'));
                        }
                        else {
                            //read dosage data
                            geno = strtod(p, NULL);
                        }
                        gsl_vector_set (genotype, ctest_idv, geno);
                        geno_mean += geno;
                        c_idv++;
                        ctest_idv++;
                    }
                }
                pch = (nch == NULL) ? NULL : nch+1;
            }
        
        geno_mean/=(double)(ni_test-n_miss);
        
        for (size_t i=0; i < ni_test; ++i) {
                if (genotype_miss[i]) {geno=geno_mean; gsl_vector_set (genotype, i, geno);}
                // do not center genotype data in UCHAR**
                else { geno = gsl_vector_get (genotype, i);}
                geno_uchar[i] = DoubleToUchar(geno);
                //UtX[ctest_snp][i] = DoubleToUchar(geno);
                //if (ctest_snp==0 && i < 10) cout << geno << ":" << (int)geno_uchar[i] << ", ";
            }
        gsl_vector_add_constant(genotype, -geno_mean); // center genotype gsl_vector here
            
            compressedBufferSize = BufferSize;
            result = compress(TempCompBuffer, &compressedBufferSize, geno_uchar, sourceBufferSize);
            if (result != Z_OK) {
                zerr(result);
                exit(-1);
            }
            else {
                UtX[ctest_snp] = (uchar*)malloc(compressedBufferSize);
                memcpy(UtX[ctest_snp], TempCompBuffer, compressedBufferSize);
                CompBuffSizeVec.push_back(compressedBufferSize);
                
                // UnCompBufferSize=sourceBufferSize;
                //  result = uncompress(TempBuffer, &UnCompBufferSize, UtX[c_snp],compressedBufferSize);
                //  if(c_snp < 10)  {
                //    zerr(result);
                //cout << "uncompressed buffer size = " << UnCompBufferSize << endl;
                //  PrintVector(TempBuffer, 10);
                // }
                // cout << "compressed Buffer size = " << compressedBufferSize << endl;
            }
            
            //JY add
            gsl_blas_ddot(genotype, genotype, &vtx);
            if(vtx < 0.00000001)
            {cout << "snp has x'x = " << setprecision(9) << vtx;}
            
            if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
            
            c_snp++;
            ctest_snp++;
        }
    }
    cout << "ctest_snp = " << c_snp << "; ns_test = " << ns_test << endl;
    
    if (calc_K==true) {
        gsl_matrix_scale (K, 1.0/(double)ns_test);
        
        for (size_t i=0; i<genotype->size; ++i) {
            for (size_t j=0; j<i; ++j) {
                geno=gsl_matrix_get (K, j, i);
                gsl_matrix_set (K, i, j, geno);
            }
        }
    }
    
    free(TempBuffer);
    free(TempCompBuffer);
    gsl_vector_free(genotype);
    delete [] geno_uchar;
    infile.clear();
    infile.close();
    
    cout << "read vcf file second time success ... \n" ;
    return true;
}

/* bool ReadFile_vcf (const string &file_vcf, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** UtX, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, vector<String> &InputSampleID, StringIntHash &sampleID2vcfInd, const string &file_sample)
{
    size_t ni_total = indicator_idv.size();
    
    vector<uint> SampleVcfPos;
    for (size_t i=0; i<ni_total; i++) {
        uint vcfpos = (uint)sampleID2vcfInd.Integer(InputSampleID[i]);
        SampleVcfPos.push_back(vcfpos);
    }
    
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    
    // Set to only store the GT genotype field.
    VcfRecordGenotype::addStoreField("GT");
    VcfRecordGenotype::addStoreField("EC");
    
    // Open the VCF file & read the header.
    
    if(!inFile.open(file_vcf.c_str(), header, file_sample.c_str(), NULL, NULL))
    {cerr<<"error reading genotype file:"<<file_vcf<<endl; exit(1);}
    
    if (calc_K==true) {gsl_matrix_set_zero (K);}
    
    gsl_vector *genotype=gsl_vector_alloc (ni_test);
    
    double geno, geno_mean;
    size_t n_miss;
    size_t c_idv=0, c_snp=0, ctest_snp = 0;
    
    while(inFile.readRecord(record))
    {
        if (!indicator_snp[c_snp]) {c_snp++; continue;}
        
        c_idv=0; geno_mean=0; n_miss=0;
        
        vector<bool> genotype_miss(ni_test, 0);
        
        for (size_t j=0; j < ni_total; ++j)
        {
            if (!indicator_idv[j]) {continue;}
            geno = getDoubleDosageFromRecord(record, SampleVcfPos[j]);
            if (geno == -9.0) {genotype_miss[c_idv]=1; n_miss++;}
            else {
                gsl_vector_set (genotype, c_idv, geno);
                geno_mean += geno;
            }
            c_idv++;
        }
        geno_mean/=(double)(ni_test-n_miss);
        
        for (size_t i=0; i < ni_test; ++i) {
            if (genotype_miss[i]) {geno=geno_mean;} // do not center genotype data in UCHAR**
            else { geno = gsl_vector_get (genotype, i);}
            
            gsl_vector_set (genotype, i, geno);
            UtX[ctest_snp][i] = DoubleToUchar(geno);
            // if (ctest_snp==0 && i < 10) cout << geno << ":" << (int)UtX[ctest_snp][i] << ", ";
        }
        gsl_vector_add_constant(genotype, -geno_mean); // center genotype gsl_vector here
        
        c_snp++;
        ctest_snp++;
        
        if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
    }
    
    if (calc_K==true) {
        gsl_matrix_scale (K, 1.0/(double)ns_test);
        
        for (size_t i=0; i<genotype->size; ++i) {
            for (size_t j=0; j<i; ++j) {
                geno=gsl_matrix_get (K, j, i);
                gsl_matrix_set (K, i, j, geno);
            }
        }
    }
    
    gsl_vector_free (genotype);
    cout << "read vcf file second time success ... \n" ;
    inFile.close();
    return true;
} */

// read multiple vcf files second time
bool ReadFile_vcfs (const string &file_vcfs, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar ** UtX, const uint ni_test, const uint ns_test, gsl_matrix *K, const bool calc_K, vector<string> &InputSampleID, StringIntHash &sampleID2vcfInd, const string &file_sample)
{
    size_t ni_total = indicator_idv.size();
    vector<uint> SampleVcfPos;
    //for (size_t i=0; i < ni_total; i++) {
       // uint vcfpos = (uint)sampleID2vcfInd.Integer(InputSampleID[i]);
      //  SampleVcfPos.push_back(vcfpos);
    //}
    
    cout << "start reading vcf files for second time ...\n";
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    // Set to only store the GT genotype field.
    VcfRecordGenotype::addStoreField("GT");
    VcfRecordGenotype::addStoreField("EC");
    ifstream infile2(file_vcfs.c_str(), ifstream::in);
    string vcf_file;
    int file_num = 0;
    
    if (calc_K==true) {gsl_matrix_set_zero (K);}
    gsl_vector *genotype=gsl_vector_alloc (ni_test);
    double geno, geno_mean;
    size_t n_miss;
    size_t c_idv=0, c_snp=0, ctest_snp = 0;
    
    while (getline(infile2, vcf_file)) {
        
        if(!inFile.open(vcf_file.c_str(), header, file_sample.c_str(), NULL, NULL))
        {cerr<<"error reading genotype file:"<<vcf_file<<endl; exit(1);}
        
        while(inFile.readRecord(record))
        {
            if (!indicator_snp[c_snp]) {c_snp++; continue;}
            
            c_idv=0; geno_mean=0; n_miss=0;
            vector<bool> genotype_miss(ni_test, 0);
            for (size_t j=0; j < ni_total; ++j)
            {
                if (!indicator_idv[j]) {continue;}
                geno = getDoubleDosageFromRecord(record, SampleVcfPos[j]);
                if (geno == -9.0) {genotype_miss[c_idv]=1; n_miss++;}
                else {
                    gsl_vector_set (genotype, c_idv, geno);
                    geno_mean += geno;
                }
                c_idv++;
            }
            geno_mean/=(double)(ni_test-n_miss);
            
            for (size_t i=0; i < ni_test; ++i) {
                if (genotype_miss[i]) {geno=geno_mean;} // do not center genotype data in UCHAR**
                else { geno = gsl_vector_get (genotype, i);}
                
                gsl_vector_set (genotype, i, geno);
                UtX[ctest_snp][i] = DoubleToUchar(geno);
            }
            gsl_vector_add_constant(genotype, -geno_mean); // center genotype gsl_vector here
            
            c_snp++;
            ctest_snp++;
            
            if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
        }
        
        inFile.close();
        file_num++;
    }
    infile2.close();
    
    if (calc_K==true) {
        gsl_matrix_scale (K, 1.0/(double)ns_test);
        
        for (size_t i=0; i<genotype->size; ++i) {
            for (size_t j=0; j<i; ++j) {
                geno=gsl_matrix_get (K, j, i);
                gsl_matrix_set (K, i, j, geno);
            }
        }
    }
    
    gsl_vector_free (genotype);
    cout << "read vcf file second time success ... \n" ;
    
    return true;
}



//Read bimbam mean genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_geno (const string &file_geno, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar **UtX, gsl_matrix *K, const bool calc_K, size_t ni_test, size_t ns_test)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
    //	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
	gsl_vector *genotype=gsl_vector_alloc (ni_test);
	gsl_vector *genotype_miss=gsl_vector_alloc (ni_test);
	double geno, geno_mean;
	size_t n_miss;
	
	uint ni_total= indicator_idv.size();
	uint ns_total= indicator_snp.size();
	
	size_t c_idv=0, c_snp=0;
	
	for (size_t i=0; i<ns_total; ++i) {
		!safeGetline(infile, line).eof();
		if (indicator_snp[i]==0) {continue;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		c_idv=0; geno_mean=0; n_miss=0;
		gsl_vector_set_zero (genotype_miss);
		for (uint j=0; j<ni_total; ++j) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[j]==0) {continue;}
			
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++;}
			else {
				geno=atof(ch_ptr);
				gsl_vector_set (genotype, c_idv, geno);
				geno_mean+=geno;
			}
			c_idv++;
		}
		
		geno_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {geno=geno_mean;}
            UtX[c_snp][i] = DoubleToUchar(geno);
			gsl_vector_set (genotype, i, (geno-geno_mean));
		}
		
		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
	}
	
	if (calc_K==true) {
		gsl_matrix_scale (K, 1.0/(double)ns_test);
		
		for (size_t i=0; i<genotype->size; ++i) {
			for (size_t j=0; j<i; ++j) {
				geno=gsl_matrix_get (K, j, i);
				gsl_matrix_set (K, i, j, geno);
			}
		}
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	
	infile.clear();
	infile.close();
	
	return true;
}


/* bool ReadFile_geno (const string &file_geno, vector<bool> &indicator_idv, vector<bool> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K)
 {
 igzstream infile (file_geno.c_str(), igzstream::in);
 //	ifstream infile (file_geno.c_str(), ifstream::in);
 if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
 
 string line;
 char *ch_ptr;
 
 if (calc_K==true) {gsl_matrix_set_zero (K);}
 
 gsl_vector *genotype=gsl_vector_alloc (UtX->size1);
 gsl_vector *genotype_miss=gsl_vector_alloc (UtX->size1);
 double geno, geno_mean;
 size_t n_miss;
 
 uint ni_total= indicator_idv.size();
 uint ns_total= indicator_snp.size();
 uint ni_test=UtX->size1;
 uint ns_test=UtX->size2;
 
 int c_idv=0, c_snp=0;
 
 for (uint i=0; i<ns_total; ++i) {
 !safeGetline(infile, line).eof();
 if (indicator_snp[i]==0) {continue;}
 
 ch_ptr=strtok ((char *)line.c_str(), " , \t");
 ch_ptr=strtok (NULL, " , \t");
 ch_ptr=strtok (NULL, " , \t");
 
 c_idv=0; geno_mean=0; n_miss=0;
 gsl_vector_set_zero (genotype_miss);
 for (uint j=0; j<ni_total; ++j) {
 ch_ptr=strtok (NULL, " , \t");
 if (indicator_idv[j]==0) {continue;}
 
 if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++;}
 else {
 geno=atof(ch_ptr);
 gsl_vector_set (genotype, c_idv, geno);
 geno_mean+=geno;
 }
 c_idv++;
 }
 
 geno_mean/=(double)(ni_test-n_miss);
 
 for (size_t i=0; i<genotype->size; ++i) {
 if (gsl_vector_get (genotype_miss, i)==1) {geno=0;}
 else {geno=gsl_vector_get (genotype, i); geno-=geno_mean;}
 
 gsl_vector_set (genotype, i, geno);
 gsl_matrix_set (UtX, i, c_snp, geno);
 }
 
 if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
 
 c_snp++;
 }
 
 if (calc_K==true) {
 gsl_matrix_scale (K, 1.0/(double)ns_test);
 
 for (size_t i=0; i<genotype->size; ++i) {
 for (size_t j=0; j<i; ++j) {
 geno=gsl_matrix_get (K, j, i);
 gsl_matrix_set (K, i, j, geno);
 }
 }
 }
 
 gsl_vector_free (genotype);
 gsl_vector_free (genotype_miss);
 
 infile.clear();
 infile.close();
 
 return true;
 } */


//Read bimbam mean genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_bed (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, uchar **UtX, gsl_matrix *K, const bool calc_K, size_t ni_test, size_t ns_test, vector <size_t> &CompBuffSizeVec)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
	
	char ch[1];
	bitset<8> b;
    double vtx;
	
	size_t ni_total= indicator_idv.size();
	size_t ns_total= indicator_snp.size();
	size_t n_bit;
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}
	
	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
    cout << "ni_test = " << ni_test << endl;
	gsl_vector *genotype = gsl_vector_alloc (ni_test);
    uchar *geno_uchar = new uchar[ni_test];
    size_t sourceBufferSize = (ni_test) * sizeof(uchar);
   // size_t UnCompBufferSize=sourceBufferSize;
    
    const size_t BufferSize = (size_t)(compressBound(sourceBufferSize));
    uchar * TempCompBuffer = (uchar*)malloc(BufferSize);
    uchar * TempBuffer = (uchar*)malloc(sourceBufferSize);

    size_t compressedBufferSize = BufferSize;
    cout << "Source Buffer Size = " << sourceBufferSize << "; Comp Buffer Bound = " << BufferSize  << endl;
    
	CompBuffSizeVec.clear();
    
	double geno, geno_mean;
	size_t n_miss;
	size_t c_idv=0, c_snp=0, c=0;
    int result;
	
	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
		if (indicator_snp[t]==0) {continue;}
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		c_idv=0; geno_mean=0.0; n_miss=0; c=0;
		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c == (size_t)ni_total) {break;}
				if (indicator_idv[c]==0) {c++; continue;}
				c++;
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); geno_mean+=2.0;}
					else {gsl_vector_set(genotype, c_idv, 1.0); geno_mean+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); geno_mean+=0.0;}
					else {gsl_vector_set(genotype, c_idv, -9.0); n_miss++;}
				}
				c_idv++;
			}
		}
        if (n_miss > 0) cout << "n_miss = " << n_miss << endl;
		if(c_idv != (size_t)ni_test) cout << "# of readed individuals not equal to ni_test \n";
        
		geno_mean/=(double)(ni_test-n_miss);
        if(geno_mean < 0.00000001) cout << "SNP_" << c_snp << "has geno_mean =" << geno_mean << endl;
        
		for (size_t i=0; i<genotype->size; ++i) {
			geno=gsl_vector_get (genotype, i);
			if (geno==-9.0) {geno=geno_mean;}
            geno_uchar[i] = DoubleToUchar(geno);
            //UtX[c_snp][i]=DoubleToUchar(geno);
            gsl_vector_set (genotype, i, (geno-geno_mean));
		}
        
        compressedBufferSize = BufferSize;
      //  if(c_snp < 10) {
        //    cout << "source geno: \n";
          //  PrintVector(geno_uchar, 10);
       // }
        result = compress(TempCompBuffer, &compressedBufferSize, geno_uchar, sourceBufferSize);
        if (result != Z_OK) {
            zerr(result);
            exit(-1);
        }
        else {
            UtX[c_snp] = (uchar*)malloc(compressedBufferSize);
            memcpy(UtX[c_snp], TempCompBuffer, compressedBufferSize);
            CompBuffSizeVec.push_back(compressedBufferSize);
            
           // UnCompBufferSize=sourceBufferSize;
          //  result = uncompress(TempBuffer, &UnCompBufferSize, UtX[c_snp],compressedBufferSize);
          //  if(c_snp < 10)  {
            //    zerr(result);
                //cout << "uncompressed buffer size = " << UnCompBufferSize << endl;
              //  PrintVector(TempBuffer, 10);
           // }
           // cout << "compressed Buffer size = " << compressedBufferSize << endl;
        }
        
        //JY add
        gsl_blas_ddot(genotype, genotype, &vtx);
        if(vtx < 0.00000001)
        {cout << "snp has x'x = " << setprecision(9) << vtx;}
		//JY add
		
		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
	}
    cout << "compressed Buffer size = " << compressedBufferSize << endl;
    cout << "CompBuffSizeVec length = " << CompBuffSizeVec.size() << endl;
    
	if(c_snp != (size_t)ns_test) cout <<"# of readed SNP not equal to ns_test \n";
	
	if (calc_K==true) {
		gsl_matrix_scale (K, 1.0/(double)ns_test);
		
		for (size_t i=0; i<genotype->size; ++i) {
			for (size_t j=0; j<i; ++j) {
				geno=gsl_matrix_get (K, j, i);
				gsl_matrix_set (K, i, j, geno);
			}
		}
	}
	
    free(TempBuffer);
    free(TempCompBuffer);
	gsl_vector_free (genotype);
    delete [] geno_uchar;
    
	infile.clear();
	infile.close();
	
	return true;
}


/*bool ReadFile_bed (const string &file_bed, vector<bool> &indicator_idv, vector<bool> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K)
 {
 ifstream infile (file_bed.c_str(), ios::binary);
 if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
 
 char ch[1];
 bitset<8> b;
 
 int ni_total=(int)indicator_idv.size();
 int ns_total=(int)indicator_snp.size();
 int ni_test=UtX->size1;
 int ns_test=UtX->size2;
 int n_bit;
 
 double vtx;
 
 if (ni_total%4==0) {n_bit=ni_total/4;}
 else {n_bit=ni_total/4+1;}
 
 
 //print the first three majic numbers
 for (int i=0; i<3; ++i) {
 infile.read(ch,1);
 b=ch[0];
 }
 
 if (calc_K==true) {gsl_matrix_set_zero (K);}
 
 gsl_vector *genotype=gsl_vector_alloc (UtX->size1);
 
 double geno, geno_mean;
 size_t n_miss;
 size_t c_idv=0, c_snp=0, c=0;
 
 //start reading snps and doing association test
 for (int t=0; t<ns_total; ++t) {
 if (indicator_snp[t]==0) {continue;}
 infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
 
 //read genotypes
 c_idv=0; geno_mean=0.0; n_miss=0; c=0;
 for (int i=0; i<n_bit; ++i) {
 infile.read(ch,1);
 b=ch[0];
 for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
 if ((i==(n_bit-1)) && c == (size_t)ni_total) {break;}
 if (indicator_idv[c]==0) {c++; continue;}
 c++;
 
 if (b[2*j]==0) {
 if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); geno_mean+=2.0;}
 else {gsl_vector_set(genotype, c_idv, 1.0); geno_mean+=1.0;}
 }
 else {
 if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); geno_mean+=0.0;}
 else {gsl_vector_set(genotype, c_idv, -9.0); n_miss++;}
 }
 c_idv++;
 }
 }
 if (n_miss > 0) cout << "n_miss = " << n_miss << endl;
 if(c_idv != (size_t)ni_test) cout << "# of readed individuals not equal to ni_test \n";
 
 geno_mean/=(double)(ni_test-n_miss);
 if(geno_mean == 0) cout << "SNP_" << c_snp << "has geno_mean = 0" << endl;
 for (size_t i=0; i<genotype->size; ++i) {
 geno=gsl_vector_get (genotype, i);
 if (geno==-9) {geno=0.0;}
 else {geno-=geno_mean;}
 
 gsl_vector_set (genotype, i, geno);
 gsl_matrix_set (UtX, i, c_snp, geno);
 }
 
 //JY add
 gsl_blas_ddot(genotype, genotype, &vtx);
 if(vtx < 0.00000001)
 {cout << "snp has x'x = " << setprecision(9) << vtx;}
 //JY add
 
 if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
 
 c_snp++;
 }
 if(c_snp != (size_t)ns_test) cout <<"# of readed SNP not equal to ns_test \n";
 
 if (calc_K==true) {
 gsl_matrix_scale (K, 1.0/(double)ns_test);
 
 for (size_t i=0; i<genotype->size; ++i) {
 for (size_t j=0; j<i; ++j) {
 geno=gsl_matrix_get (K, j, i);
 gsl_matrix_set (K, i, j, geno);
 }
 }
 }
 
 gsl_vector_free (genotype);
 infile.clear();
 infile.close();
 
 return true;
 } */





bool ReadFile_est (const string &file_est, const vector<size_t> &est_column, map<string, double> &mapRS2est)
{
	mapRS2est.clear();
	
	ifstream infile (file_est.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening estimated parameter file: "<<file_est<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	double alpha, beta, gamma, d;
	
	//header
	getline(infile, line);
	
	size_t n=*max_element(est_column.begin(), est_column.end());
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		
		alpha=0.0; beta=0.0; gamma=1.0;
		for (size_t i=0; i<n+1; ++i) {
			if (i==est_column[0]-1) {rs=ch_ptr;}
			if (i==est_column[1]-1) {alpha=atof(ch_ptr);}
			if (i==est_column[2]-1) {beta=atof(ch_ptr);}
			if (i==est_column[3]-1) {gamma=atof(ch_ptr);}
			if (i<n) {ch_ptr=strtok (NULL, " \t");}
		}
		
		d=alpha+beta*gamma;
		
		if (mapRS2est.count(rs)==0) {
			mapRS2est[rs]=d;
		}
		else {
			cout<<"the same SNP occurs more than once in estimated parameter file: "<<rs<<endl; return false;
		}
	}
	
	infile.clear();
	infile.close();
	return true;
}



bool CountFileLines (const string &file_input, size_t &n_lines)
{
	igzstream infile (file_input.c_str(), igzstream::in);
	//ifstream infile (file_input.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open file: "<<file_input<<endl; return false;}
    
	n_lines=count(istreambuf_iterator<char>(infile), istreambuf_iterator<char>(), '\n');
	infile.seekg (0, ios::beg);
	
	return true;
}



//Read gene expression file
bool ReadFile_gene (const string &file_gene, vector<double> &vec_read, vector<SNPINFO> &snpInfo, size_t &ng_total)
{
	vec_read.clear();
	ng_total=0;
	
	ifstream infile (file_gene.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open gene expression file: "<<file_gene<<endl; return false;}
	
	string line;
	char *ch_ptr;
	string rs;
	
	size_t n_idv=0, t=0;
    vector<bool> indicator_func_temp;
    vector<double> weight_temp;
    
	//header
	getline(infile, line);
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		
		ch_ptr=strtok (NULL, " , \t");
		
		t=0;
		while (ch_ptr!=NULL) {
			if (ng_total==0) {
				vec_read.push_back(0);
				t++;
				n_idv++;
			} else {
				vec_read[t]+=atof(ch_ptr);
				t++;
			}
			
			ch_ptr=strtok (NULL, " , \t");
		}
		
		if (t!=n_idv) {cout<<"error! number of columns doesn't match in row: "<<ng_total<<endl; return false;}
		
		SNPINFO sInfo={"-9", rs, -9, -9, "-9", "-9", -9, -9, -9, indicator_func_temp, weight_temp, 0.0};
		snpInfo.push_back(sInfo);
		
		ng_total++;
	}
	
    infile.clear();
	infile.close();
	
	return true;
}

