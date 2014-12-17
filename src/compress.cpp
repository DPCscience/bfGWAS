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








