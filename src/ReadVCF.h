#include "VcfFileReader.h"
#include "StringBasics.h"
#include "StringHash.h"
#include "MemoryAllocators.h"

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#include "/net/fantasia/home/yjingj/software/zlib-1.2.8/zlib.h"
#include "/net/fantasia/home/yjingj/software/gzstream/gzstream.h"
//#include "/net/fantasia/home/yjingj/software/miniz/miniz.c"

#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include <fstream>
//#include <string>
#include <sstream>

#include "param.h"

typedef short int int16;
typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;

//#define my_max(a,b) (((a) > (b)) ? (a) : (b))
//#define my_min(a,b) (((a) < (b)) ? (a) : (b))

//#define BUF_SIZE (1024 * 1024)
//static uint8 s_inbuf[BUF_SIZE];
//static uint8 s_outbuf[BUF_SIZE];




struct genotypeMatrix
{
    std::vector<String> sampleIDs;
    std::vector<genMarker> markers;
    
    uchar** genotypes;
    uint genotypesSize;
    
    void initializeMatrix(const char* filename);
    
};


bool print(const char* description, genotypeMatrix& info, uint numSamples, uint numMarkers);

bool print(const char* description, uchar **genotypes, uint numMarkers, uint numSamples, vector<String> &sampleIDs);


float StringToFloat(const char* s);
float UcharToFloat(const uchar c);
double UcharToDouble(const uchar c);

uchar FloatToUchar(const float doseage);
uchar IntToUchar(const int intc);
uchar DoubleToUchar(const double doseage);

void getGTgslVec(uchar ** X, gsl_vector *xvec, size_t marker_i, const size_t ni_test, const size_t ns_test);
bool getGTgslMat(uchar ** X, gsl_vector *Xgsl, vector<size_t> marker_idx, const size_t ni_test, const size_t ns_test);

uchar getUcharDosageFromRecord(VcfRecord &record, const uint smNum);

