#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "zlib.h"
#include "zconf.h"

//#include "stdafx.h"   // not actually needed
//#define ZLIB_WINAPI   // actually actually needed (for linkage)

//#include "windows.h"  // get BYTE et al.
//#pragma comment(lib, "zlibwapi.lib") // for access to the DLL


#ifndef uchar
typedef unsigned char uchar;
#endif
typedef unsigned char BYTE;



void zerr(int ret);

int CompressGenoVec(uchar * sourceBuffer, size_t sourceBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize);

int DecompressGenoVec(uchar * UnCompBuffer, size_t &UnCompBufferSize, uchar * compressedBuffer, size_t compressedBufferSize);

int GetMaxCompressedLen(int nLenSrc );

int CompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );

int UncompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );
