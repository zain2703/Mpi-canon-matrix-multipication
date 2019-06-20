#ifndef __UTILS_H__
#define __UTILS_H__

#include "mpi.h"

void* g_malloc(size_t size);

void* g_calloc( size_t num, size_t size);

void* g_realloc( void *ptr, size_t size);

void readMatrices( const char* fileName, double** A, double** B, int* N  );

void writeMatrix( const char* fileName, double* numbers, int N);

void exitError( const char* msg, ... );

#ifdef DEBUG
    void trace( const char* msg, ... );
    void traceVector( double* numbers, int backetSize );
#else
    #define trace(expr, ...)                ((void)0)
    #define traceVector( expr1, expr2 )     ((void)0)
#endif

#endif
