#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "zlib.h"

#ifndef uchar
typedef unsigned char uchar;
#endif

void zerr(int ret);

int CompressGenoVec(uchar * sourceBuffer, size_t &sourceBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize);

int DecompressGenoVec(uchar * UnCompBuffer, size_t &UnCompBufferSize, uchar * compressedBuffer, size_t &compressedBufferSize);
