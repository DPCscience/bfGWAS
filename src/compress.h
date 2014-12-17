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

int GetMaxCompressedLen(int nLenSrc );

int CompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );

int UncompressData( const BYTE* abSrc, int nLenSrc, BYTE* abDst, int nLenDst );
