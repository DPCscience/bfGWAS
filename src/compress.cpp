#include "compress.h"


//size_t sourceBufferSize = ni_test * sizeof(uchar);
//uchar* sourceBuffer = (uchar *)malloc(sourceBufferSize);
//size_t compressedBufferSize = compressBound(sourceBufferSize);
/* report a zlib or i/o error */

void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
        case Z_ERRNO:
            if (ferror(stdin))
                fputs("error reading stdin\n", stderr);
            if (ferror(stdout))
                fputs("error writing stdout\n", stderr);
            break;
        case Z_STREAM_ERROR:
            fputs("invalid compression level\n", stderr);
            break;
        case Z_DATA_ERROR:
            fputs("invalid or incomplete deflate data\n", stderr);
            break;
        case Z_MEM_ERROR:
            fputs("out of memory\n", stderr);
            break;
        case Z_VERSION_ERROR:
            fputs("zlib version mismatch!\n", stderr);
    }
}

int CompressGenoVec(uchar * sourceBuffer, size_t &sourceBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize){
    
    compressedBuffer = (uchar *)malloc(compressedBufferSize);
    int result = compress(compressedBuffer, &compressedBufferSize, sourceBuffer, sourceBufferSize);
    
    if (result != Z_OK) {
        zerr(result);
        return result;
    }
    else return 1;
    
    uchar * newBuffer = (uchar *)realloc((void*)compressedBuffer, compressedBufferSize);
    if (newBuffer != NULL) {
        compressedBuffer = newBuffer;
        return 1;
    }
    else {
        std::cerr << "Error compressing buffer\n";
        exit(-1);
    }
    
}


int DecompressGenoVec(uchar * UnCompBuffer, size_t &UnCompBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize){
    
    UnCompBuffer = (uchar *)malloc(UnCompBufferSize);
    int result = uncompress(UnCompBuffer, &UnCompBufferSize, compressedBuffer, compressedBufferSize);

    if (result != Z_OK) {
        zerr(result);
        return result;
    }
    else return 1;
    
}

/*
bool CompressGenotype(uchar ** X, size_t ns_test, size_t ni_test, std::vector<uchar*> Xcompressed){

    size_t sourceBuffSize = ni_test * sizeof(uchar);
    
    size_t compressedBuffSize = compressBound(sourceBuffSize);
    
    compressedBuffer = (uchar *)malloc(comprssedBufferSize);
    int result = compress(compressedBuffer, &compressedBufferSize, sourceBuffer, sourceBufferSize);
    
    if (result != Z_OK) {
        zerr(result);
        return result;
    }
    else return 1;
    
    uchar * newBuffer = (uchar *)realloc((void*)compressedBuffer, compressedBufferSize);
    if (newBuffer != NULL) {
        compressedBuffer = newBuffer;
        return 1;
    }
    else {
        std::cerr << "Error compressing buffer\n";
        exit(-1);
    }

} */










