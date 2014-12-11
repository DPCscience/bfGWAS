#include "ReadVCF.h"



void genotypeMatrix::initializeMatrix(const char* filename)
{
    genotypesSize = 0;
    
    IFILE inputFile = ifopen(filename, "r");
    if(inputFile == NULL)
    {
        std::cerr << "Unable to open " << filename << "\n";
        exit(1);
    }
    
    while(inputFile->ifgetc() == '#')
    {

        if(inputFile->discardLine() < 0)
        {
            break;
        }
    }
    
    // Count data lines
    while(inputFile->discardLine() >= 0 && genotypesSize < 999)
    {
        ++genotypesSize;
    }
    // Increase genotypesSize size by 1 in case there
    // is not a new line at end of last marker
    ++genotypesSize;
    std::cout << "Total # of markers to load = " << genotypesSize << "\n";
    ifclose(inputFile);
    
    ////////////////////////////////////////
    // Now we have the number of markers.
    // Pass 2 - store marker information.
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    
    // Set to only store the GT genotype field.
    VcfRecordGenotype::addStoreField("GT");
    VcfRecordGenotype::addStoreField("EC");
    
    // Open the VCF file & read the header.
    inFile.open(filename, header);
    
    uint numSamples = (uint)header.getNumSamples();
    if (numSamples>100) {
        numSamples=100;
    }
    std::cout << "Total # of Samples to load = " << numSamples << "\n";
    // Store the sample names
    for(uint i = 0; i < numSamples; i++)
    {
        sampleIDs.push_back(header.getSampleName(i));
    }
    
    ////////////
    // Allocate the genotypes Matrix.
    //Testing a subset data
    genotypes = AllocateUCharMatrix(genotypesSize, numSamples);
    InitialUCharMatrix(genotypes, genotypesSize, numSamples, UCHAR_MAX);
    
    // Read each record
    uint markerIndex = 0;
    genMarker temp_genMarker;
    
    while(inFile.readRecord(record) && markerIndex < genotypesSize)
    {
        // Add the marker
        temp_genMarker.iniRecord(record);
        markers.push_back(temp_genMarker);
        // Add each individual's genotype
        
        for(uint smNum = 0; smNum < numSamples; smNum++)
        {
            genotypes[markerIndex][smNum] = getUcharDosageFromRecord(record, smNum);
        }
        //std::cout << "\n";
        ++markerIndex;
    }
    
}

double getDoubleDosageFromRecord(VcfRecord& record, const uint smNum)
{
    //Read EC, dosage data from vcf in hex string
    double dosage = -9.0;
    
    const std::string* ecStrPtr = record.getGenotypeInfo().getString("EC", smNum);
    if(ecStrPtr != NULL)
    {
        //cout << *ecStrPtr << " ";
        dosage = strtod(ecStrPtr->c_str(), NULL);
    }
    
    //read GT
    else{
        int alleleIndex;
        dosage = 0.0;
        for(int i = 0; i < 2; i++)
        {
            alleleIndex = record.getGT(smNum, i);
            if(alleleIndex == VcfGenotypeSample::MISSING_GT)
            { dosage = -9.0; break;}
            else if(alleleIndex > 0)
            {
                dosage += 1.0;
            }
        }
    }

    if (dosage < 0.0 || dosage > 2.0) {
        dosage = -9.0;
    }
    
    return dosage;
    
}


uchar getUcharDosageFromRecord(VcfRecord& record, const uint smNum)
{
    //Read EC, dosage data from vcf in hex string
    uchar c;
    int intc;
    float dosage;
    
    const std::string* ecStrPtr = record.getGenotypeInfo().getString("EC", smNum);
    if(ecStrPtr != NULL)
    {
        dosage = strtof(ecStrPtr->c_str(), NULL);
        c=FloatToUchar(dosage);
    }

//read GT
else{
    intc = 0;
    int alleleIndex;
    for(int i = 0; i < 2; i++)
    {
        alleleIndex = record.getGT(smNum, i);
        if(alleleIndex == VcfGenotypeSample::MISSING_GT)
        { intc = -2; break;}
        else if(alleleIndex > 0)
        { intc += 1; }
    }
    c = IntToUchar(intc);
    }
    return c;
}

double StringToDouble(const char* s) {
    double f = strtod(s, NULL);
    return f;
}

float StringToFloat(const char* s) {
    float f = strtof(s, NULL);
    return f;
}

float UcharToFloat(const uchar c){
    if (c != UCHAR_MAX) {
        int intc = c;
        return (((float)intc) * 0.01);
    }
    else return -9.0;
}

double UcharToDouble(const uchar c){
    if (c != UCHAR_MAX) {
        int intc = c;
        return (((double)intc) * 0.01);
    }
    else return -9.0;
}

uchar FloatToUchar(const float dosage){
    if (dosage >= 0.0  && dosage <= 2.0) {
        return  ((int)(dosage*100.0));
    }
    else return UCHAR_MAX;
}

uchar DoubleToUchar(const double dosage){
    if (dosage >= 0.0  && dosage <= 2.0) {
        return ((int)(dosage*100.0));
    }
    else return UCHAR_MAX;
}

uchar IntToUchar(const int intc){
    if (intc >= 0  && intc <= 2) {
        return (intc * 100);
    }
    else return UCHAR_MAX;
}

// get genotype vector gor given column
void getGTgslVec(uchar ** X, gsl_vector *xvec, size_t marker_i, const size_t ni_test, const size_t ns_test){

    double geno, geno_mean = 0.0;
    if (marker_i < ns_test ) {
        for (size_t j=0; j<ni_test; j++) {
            geno = UcharToDouble(X[marker_i][j]);
            //cout << geno << ", ";
            if (geno < 0.0 || geno > 2.0) {
                cerr << "wrong genotype value = " << geno << endl;
                exit(1);
            }
            gsl_vector_set(xvec, j, geno);
            geno_mean += geno;
        }
        geno_mean /= (double)ni_test;
        geno_mean = -geno_mean;
        gsl_vector_add_constant(xvec, geno_mean); // center genotypes here
    }
    else {
        std::cerr << "Error return genotype vector...\n";
        exit(1);
    }
    //cout << endl;
}

//get genotype matrix for given column vector
bool getGTgslMat(uchar ** X, gsl_matrix *Xgsl, vector<size_t> marker_idx, const size_t ni_test, const size_t ns_test){
    
    size_t marker_i;
    gsl_vector *xvec = gsl_vector_alloc(ni_test);
    
    for (size_t i=0; i < marker_idx.size(); i++) {
        marker_i = marker_idx[i];
        getGTgslVec(X, xvec, marker_i, ni_test, ns_test);
        gsl_matrix_set_col(Xgsl, i, xvec);
    }
    
    gsl_vector_free(xvec);
    return 1;
}

// print sub genotype uchar array
bool print(const char* description, uchar **genotypes, uint numMarkers, uint numSamples, vector<String> &sampleIDs)
{
    std::cout << "\n\t\t\t" << description << "\n\t\t";
    std::cout << "Sample ID : ";
    std::cout << "\n";
    for(uint i = 0; i < numSamples; i++)
    {
        std::cout << "\t" << sampleIDs[i];
    }
    std::cout << "\n";
    
    uchar c;
    float fc;
    std::cout << "\n" << "Dosage GT: \n";
    for(uint i = 0; i < numMarkers; i++)
    {
      for(uint j = 0; j < numSamples; j++)
        {
            std::cout << "\t";
            if(genotypes[i][j] != UCHAR_MAX)
            {
                c = genotypes[i][j];
                fc = (float)c;
                std::cout << fc * 0.01;
            }
            else
                std::cout << "GT missing";
        }
        std::cout << "\n";
    }
    
    std::cout << "Succesfully output sub-genotype\n";
    return 1;
}

bool print(uchar **genotypes, uint numMarkers, uint numSamples)
{
    int intc;
    float fc;
    std::cout << "\n" << "Dosage GT: \n";
    for(uint i = 0; i < numMarkers; i++)
    {
        for(uint j = 0; j < numSamples; j++)
        {
            std::cout << "\t";
            intc = (int)genotypes[i][j];
            //std::cout << intc << ":" ;
            
            if(genotypes[i][j] != UCHAR_MAX)
            {
                fc = (float)(intc);
                std::cout << fc * 0.01;
            }
            else
                std::cout << "GT missing";
        }
        std::cout << "\n";
    }
    
    std::cout << "Succesfully output sub-genotype\n";
    return 1;
}


bool print(const char* description, genotypeMatrix& info, uint numMarkers, uint numSamples)
{
    std::cout << "\n\t\t\t" << description << "\n\t\t";
    std::cout << "Sample ID : ";
    std::cout << "\n";
    for(uint i = 0; i < numSamples; i++)
    {
        std::cout << "\t" << info.sampleIDs[i];
    }
    std::cout << "\n";
    
    uchar c;
    std::cout << "\n" << "Dosage GT: \n";
    for(uint i = 0; i < numMarkers; i++)
    {
        info.markers[i].printMarker();
        for(uint j = 0; j < numSamples; j++)
        {
            std::cout << "\t";
            
            if(info.genotypes[i][j] != UCHAR_MAX)
            {c = info.genotypes[i][j];
                std::cout << (float)c * 0.01;}
            else
                std::cout << "GT missing";
        }
        std::cout << "\n";
    }
    
    std::cout << "Succesfully output sub-genotype\n";
    return 1;
}


