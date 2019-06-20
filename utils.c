
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>
#include <stdarg.h>
#include <math.h>

#include <mpi.h>

#include "utils.h"

void exitError( const char* msg, ... )
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    char buf[255] = "";
    sprintf( buf, "[%d] ERROR: %s", rank, msg );

    va_list args;
    va_start( args, msg );
        printf("\n");
        vprintf( buf, args ); 
        printf("\n\n");
    va_end( args );

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE );
} // exitError


#ifdef DEBUG
void trace( const char* msg, ... )
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    char buf[255] = "";
    sprintf( buf, "[%d] %s", rank, msg );

    va_list args;
    va_start( args, msg );
        vprintf( buf, args ); 
    va_end( args );
} // trace

void traceVector( double* numbers, int backetSize )
{
    for( int i = 0; i < backetSize; ++i )
        trace( "numbers[%d] = %0.6f\n", i, numbers[i] );
} // traceVector
#endif


void* g_malloc(size_t size)
{
    void* p = malloc( size );
    if( p == NULL )
        exitError( "(malloc:%d): Now you've done it. Out of memory!", size );

    return( p );
} // g_malloc


void* g_calloc( size_t num, size_t size)
{
    void* p = calloc( num, size );
    if( p == NULL )
        exitError( "(calloc:%d): Now you've done it. Out of memory!", size );

    return( p );
} // g_calloc


void* g_realloc( void *ptr, size_t size)
{
    // Avoid realloc per implementation defined behavior when size is zero.
    if( size == 0 )
        return( ptr );

    void* p = realloc( ptr, size );
    if( p == NULL )
        exitError( "(realloc:%d): Now you've done it. Out of memory!", size );

    return( p );
} // g_realloc

 void readMatrices( const char* fileName, double** A, double** B, int* N )
{
    FILE* fp = fopen( fileName, "r" ); 
    if( fp == NULL ) 
        return;
    
    char buf[100] = "";
    double* numbers = NULL;
    int  cc = 0;
    // Read number of integers.
    if( fscanf( fp, "%s\n\n", buf ) == 1 )
    {
		// Read Matrix Dimension.
        cc = atoi( buf ); 
        *N = cc; 		
		numbers = (double* )g_calloc( (*N) * (*N), sizeof(double)); 
		// Read A
        for(int i=0; i<(*N)*(*N); i++) 
		{ fscanf(fp, "%lf", &(numbers[i]));}
        *A = numbers;
        // Read B
        numbers = (double* )g_calloc( (*N) * (*N), sizeof(double) );
        for(int i=0; i<(*N)*(*N); i++) 
		{fscanf(fp, "%lf", &(numbers[i]));}
        *B = numbers;
    }
    fclose( fp );
} // readMatrices


void writeMatrix( const char* fileName, double* numbers, int N)
{
    FILE* fp = fopen(fileName, "w+" ); 
    if( fp == NULL ) 
        return;
    for( int i = 0; i < N; ++i )
    {
        for( int j = 0; j < N; ++j )
            fprintf( fp, "%lf ", numbers[i*N +j] );
        fprintf( fp, "\n" );
    }
    fclose( fp );
} // writeMatrix