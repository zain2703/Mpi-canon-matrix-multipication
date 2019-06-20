#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "utils.h"

// File parameters.
char* inputFile   = NULL;
char* outputFile  = NULL;

int globalRank = 0;
int globalNump = 0;

double executionTime = 0.0;

MPI_Status status;
MPI_Comm comm;

//Memory load
double *P_A, *P_B, *P_C;

void startUp( int argc, char *argv[] );
void done();
void check (int N, int P_Sroot);
void load(int N, int P_Sroot,int P_size); 
int get_rank(int row, int col, int P_Sroot); 
void shuffle(double** A, double *matrix, int n, int P_Sroot); 
void IntialShift(int P_Sroot, int P_size);
void MatrixMul(int P_part);
void Canon(int P_Sroot, int P_size, int P_part);
void Gather(double **C,int N,int P_size,int P_Sroot,int P_part);
	
int main(int argc, char *argv[])
{ 
    double *A,*B,*C;
    int N, P_Sroot,P_part,P_size;
    
	startUp(argc, argv);
	P_Sroot=sqrt(globalNump);			
	if(globalRank==0)
	{	
		readMatrices( inputFile,&A,&B,&N);
	    check(N,P_Sroot);
		shuffle(&A,A,N, P_Sroot);
	    shuffle(&B,B,N, P_Sroot);
	}
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD); 
     
	 P_part = N/P_Sroot;
	 P_size = pow(P_part,2);  

	 //load dynamic memory
	 load(N,P_Sroot,P_size); 

	executionTime = MPI_Wtime();
	
	MPI_Scatter(A,P_size,MPI_DOUBLE,P_A,P_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(B,P_size,MPI_DOUBLE,P_B,P_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	//paramters for Mpi_cartesian plan
	 int dim[2], period[2], reorder;
          dim[0]=P_Sroot; dim[1]=P_Sroot;  period[0]=1;period[1]=1;reorder=1;	  
	 MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
	 
	 //Canon Algorithem
	 IntialShift(P_Sroot,P_size);
	 Canon(P_Sroot, P_size, P_part);
	 
	//Gather from different processors
	 Gather(&C,N,P_size,P_Sroot,P_part);
	   
	// Measure time. Excluding file write. 
     executionTime = MPI_Wtime() - executionTime;
	//free dynamic memory
    	free(P_A);free(P_B);free(P_C);
		
	     if(globalRank==0)
	   {
		 	writeMatrix( outputFile, C, N); 		
	   } 
	   
	 done();
   return 0;
} // main

void startUp( int argc, char *argv[] )
{
    // Initializing MPI.
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &globalNump );
    MPI_Comm_rank( MPI_COMM_WORLD, &globalRank );

    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name( name, &len );
    trace( "Rank %d running on %s\n", globalRank, name );

    if( argc < 3 )
        exitError("Usage: matmul inputFile outputFile" );

    inputFile  = strdup( argv[1] );
    outputFile = strdup( argv[2] );   
} // startUp


void done()
{
    // Print out time for each process.
    trace("Elapsed time       : %lf ws\n", executionTime);

    // Wait for all to get here and write biggest execution time.
    double maxTime = 0.0;
    MPI_Reduce( &executionTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

    if( globalRank == 0 )
	{printf("%lf\n", maxTime);}
	MPI_Finalize();
} // done


void check (int N, int P_Sroot)
{
	if(globalNump != pow(P_Sroot,2)) {
		exitError("Error: process number is not a square");
    }
	
	if(N%P_Sroot != 0) {
        exitError("Error: n is not divisible");
    }
}

 void load(int N, int P_Sroot,int P_size) 
 {	      
	 P_A=(double *)g_malloc(P_size*sizeof(double));
	 P_B=(double *)g_malloc(P_size*sizeof(double));
	 P_C=(double *)g_calloc(P_size,sizeof(double));
	 memset( P_C, 0, P_size*sizeof(double) );
 }
 
 int get_rank(int row, int col, int P_Sroot) {
    int rank = ((row+P_Sroot)%P_Sroot)*P_Sroot+(col+P_Sroot)%P_Sroot;
    return rank;
}
 
 void shuffle(double** A, double *matrix, int n, int P_Sroot) {
    int P_part = n/P_Sroot;
	  int c=0;
	//matrix P_size array for shuffle
    double *temp = (double *)malloc((n*n)*sizeof(double));
     for(int i=0; i<P_Sroot; i++) {
        for(int j=0; j<P_Sroot; j++) {
            int ii = i*P_part*n;
            int jj = j*P_part;
            for(int k=0; k<P_part; k++) {
                for(int l=0; l<P_part; l++) {
                    temp[c++] = matrix[ii+jj+l];
                }
                ii += n;
            }
		}
    }
	*A=temp;
}
 
 void IntialShift(int P_Sroot, int P_size)
{ 
     int s, r,coord[2];
	 MPI_Cart_coords(comm, globalRank, 2, coord);	 
	for(int i=0; i<coord[0]; i++) {    
		s =get_rank(coord[0], coord[1]-1, P_Sroot);
		r  =get_rank(coord[0], coord[1]+1, P_Sroot);
		MPI_Sendrecv_replace(P_A, P_size, MPI_DOUBLE, s, 123, r, 123, MPI_COMM_WORLD, &status);
    }

    for(int j=0; j<coord[1]; j++) {
		s =get_rank(coord[0]-1, coord[1], P_Sroot);
		r  =get_rank(coord[0]+1, coord[1], P_Sroot);   
        MPI_Sendrecv_replace(P_B, P_size, MPI_DOUBLE, s, 111, r, 111, MPI_COMM_WORLD, &status); 
    }	
}


void MatrixMul(int P_part) {
    double* temp = (double *)calloc(P_part*P_part, sizeof(double));
    for(int i=0; i<P_part; i++){
        for(int j=0; j<P_part; j++){
            for(int k=0; k<P_part; k++){
                temp[i*P_part+j] += P_A[i*P_part+k]*P_B[k*P_part+j];
            }
            P_C[i*P_part+j] += temp[i*P_part+j];
        }
    }
    free(temp);
}

 void Canon(int P_Sroot, int P_size, int P_part)
{ 
     int s, r,coord[2];
	for(int i=0; i<P_Sroot; i++) {
        MatrixMul(P_part);
		MPI_Cart_coords(comm, globalRank, 2, coord);
		s =get_rank(coord[0], coord[1]-1, P_Sroot);
		r  =get_rank(coord[0], coord[1]+1, P_Sroot);
		MPI_Sendrecv_replace(P_A, P_size, MPI_DOUBLE, s, 222, r, 222, MPI_COMM_WORLD, &status);
 
		s =get_rank(coord[0]-1, coord[1], P_Sroot);
		r  =get_rank(coord[0]+1, coord[1], P_Sroot);   
        MPI_Sendrecv_replace(P_B, P_size, MPI_DOUBLE, s, 333, r, 333, MPI_COMM_WORLD, &status); 
    }	
}

void Gather(double **C,int N,int P_size,int P_Sroot,int P_part)
{
	double* matrix = (double *)malloc(N*N*sizeof(double));
	if(globalRank != 0) {
        MPI_Send(P_C, P_size, MPI_DOUBLE, 0, 555, MPI_COMM_WORLD);
       } 
	else {
		
        double* temp = (double *)calloc(P_size, sizeof(double));
        for(int i=0; i<P_Sroot; i++) {
            for(int j=0; j<P_Sroot; j++) {
                if(i==0 && j==0) {
                    memcpy(temp, P_C, P_size*sizeof(double));
                } else {
                    MPI_Recv(temp, P_size, MPI_DOUBLE, i*P_Sroot+j, 555, MPI_COMM_WORLD, &status);
                }
                int ii = i*P_part*N;
                int jj = j*P_part;
                for(int k=0; k<P_part; k++) {
                    for(int l=0; l<P_part; l++) {
                        matrix[ii+jj+l] = temp[k*P_part+l];
                    }
                    ii += N;
                }
            }
        }
        free(temp);
    } 		
	*C=matrix;
}
