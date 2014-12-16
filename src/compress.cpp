#include "compress.h"
#include "bslmm.h"

//size_t sourceBufferSize = ni_test * sizeof(uchar);
//uchar* sourceBuffer = (uchar *)malloc(sourceBufferSize);
//size_t compressedBufferSize = compressBound(sourceBufferSize);
/* report a zlib or i/o error */

void zerr(int ret)
{
   // fputs("zpipe: ", stderr);
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


/////////////////////////////

int GetMaxCompressedLen( int nLenSrc )
{
    size_t n16kBlocks = (nLenSrc+16383) / 16384; // round up any fraction of a block
    return ( nLenSrc + 6 + (n16kBlocks*5) );
}

int CompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst )
{
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= (BYTE*)abSrc;
    zInfo.next_out= abDst;
    
    int nErr, nRet= -1;
    nErr= deflateInit( &zInfo, Z_DEFAULT_COMPRESSION ); // zlib function
    if ( nErr == Z_OK ) {
        nErr= deflate( &zInfo, Z_FINISH );              // zlib function
        if ( nErr == Z_STREAM_END ) {
            nRet= zInfo.total_out;
        }
    }
    deflateEnd( &zInfo );    // zlib function
    return( nRet );
}

int UncompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst )
{
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= (BYTE*)abSrc;
    zInfo.next_out= abDst;
    
    int nErr, nRet= -1;
    nErr= inflateInit( &zInfo );               // zlib function
   // std::cout << "nErr of inflateInit =" << nErr << "; ";

    if ( nErr == Z_OK ) {
        nErr= inflate( &zInfo, Z_FINISH );     // zlib function
      //  std::cout << "nErr of inflate =" << nErr << "; ";

        if ( nErr == Z_STREAM_END ) {
            nRet= zInfo.total_out;
        }
      //  std::cout << "final nErr = " << nErr << "\n ";
    }
    inflateEnd( &zInfo );   // zlib function
    return( nRet ); // -1 or len of output
}

///////////////////////

int CompressGenoVec(uchar * sourceBuffer, size_t sourceBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize){
    
    //std::cout << "create tempBuffer with size " << compressedBufferSize << "\n";
    uchar * tempBuffer = (uchar *)malloc(compressedBufferSize);
  //  uchar * tempUncompBuffer = (uchar *)malloc(sourceBufferSize);
    
    //int tempSize = (int)compressedBufferSize;
    std::cout << "source Buffer : \n" ;
  //  PrintVector(sourceBuffer, 10);
    int ret = CompressData( sourceBuffer, (int)sourceBufferSize, tempBuffer, (int)compressedBufferSize );
    if (ret>0) {
        
        compressedBufferSize=(size_t)ret;
        compressedBuffer = (uchar *)malloc(compressedBufferSize);
        std::cout << "copy tempBuffer to compressedBuffer with size " << compressedBufferSize << "\n";
        memcpy(compressedBuffer, tempBuffer, compressedBufferSize);
        
    }
    else {std::cerr << "Error compressing buffer...\n";
        exit(-1);}
    
    //compress(tempBuffer, &compressedBufferSize, sourceBuffer, sourceBufferSize);
   // std::cout << "compress success ...\n";
   // if (result != Z_OK) {
     //   zerr(result);
       // return result;
    //}
    //else return 1;
    
    free(tempBuffer);
  //  free(tempUncompBuffer);
    return 1;
    
}


int DecompressGenoVec(uchar * UnCompBuffer, size_t &UnCompBufferSize, uchar * compressedBuffer, size_t compressedBufferSize){
    
    std::cout << "create tempBuffer with size " << UnCompBufferSize << "\n";
    
    uchar * tempBuffer = new uchar [UnCompBufferSize];
    int result = UncompressData( compressedBuffer, (int)compressedBufferSize, tempBuffer, (int) UnCompBufferSize );
    //uncompress(tempBuffer, &UnCompBufferSize, compressedBuffer, compressedBufferSize);
    std::cout << "result = " << result;
    if (result > 0) {
        UnCompBufferSize=(size_t)result;
    }
    else {
        std::cerr << "Error uncompressing buffer...\n";
        exit(-1);
    }
    std::cout << "ending uncompress with UnCompBufferSize ="<< UnCompBufferSize<< "\n";
    
   // if (result != Z_OK) {
     //   zerr(result);
       // return result;
    //}
    
    UnCompBuffer = new uchar [UnCompBufferSize];
    memcpy(UnCompBuffer, tempBuffer, UnCompBufferSize);
    delete [] tempBuffer;
    
    if (UnCompBuffer != NULL) {
        //UnCompBuffer = newBuffer;
        return 1;
    }
    else {
        std::cerr << "Error uncompressing buffer...\n";
        exit(-1);
    }
    
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










