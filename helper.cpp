/* 
 * Utilities for the Aliev-Panfilov code
 * Scott B. Baden, UCSD
 * Nov 2, 2015
 */

#include <iostream>
#include <assert.h>
// Needed for memalign
#include <malloc.h>
#include <mutex>

using namespace std;

// Barrier Class
/*
class barrier {
   int count, _NT;
   mutex arrival;
   mutex departure;

  public:
   barrier(int NT){
      arrival.unlock();
      departure.lock();
      count = 0;
      _NT = NT;
   }

   void bsync(){
      arrival.lock();
      if(++count<_NT)
        arrival.unlock();
      else
        departure.unlock();
      departure.lock();
      if(--count > 0)
        departure.unlock();
      else
        arrival.unlock();
   }
};
*/
void printMat(const char mesg[], double *E, int m, int n);

//
// Initialization
//
// We set the right half-plane of E_prev to 1.0, the left half plane to 0
// We set the botthom half-plane of R to 1.0, the top half plane to 0
// These coordinates are in world (global) coordinate and must
// be mapped to appropriate local indices when parallelizing the code
//
void init (double *E,double *E_prev,double *R,int m,int n){
    int i;

    for (i=0; i < (m+2)*(n+2); i++)
        E_prev[i] = R[i] = 0;

    for (i = (n+2); i < (m+1)*(n+2); i++) {
	int colIndex = i % (n+2);		// gives the base index (first row's) of the current index

        // Need to compute (n+1)/2 rather than n/2 to work with odd numbers
	if(colIndex == 0 || colIndex == (n+1) || colIndex < ((n+1)/2+1))
	    continue;

        E_prev[i] = 1.0;
    }

    for (i = 0; i < (m+2)*(n+2); i++) {
	int rowIndex = i / (n+2);		// gives the current row number in 2D array representation
	int colIndex = i % (n+2);		// gives the base index (first row's) of the current index

        // Need to compute (m+1)/2 rather than m/2 to work with odd numbers
	if(colIndex == 0 || colIndex == (n+1) || rowIndex < ((m+1)/2+1))
	    continue;

        R[i] = 1.0;
    }
    // We only print the meshes if they are small enough
    printMat("E_prev",E_prev,m,n);
    printMat("R",R,m,n);
}

double *alloc1D(int m,int n){
    int nx=n, ny=m;
    double *E;
    // Ensures that allocatdd memory is aligned on a 16 byte boundary
    assert(E= (double*) memalign(16, sizeof(double)*nx*ny) );
    return(E);
}

void printMat(const char mesg[], double *E, int m, int n){
    int i;
    if (m>8)
      return;
    printf("%s\n",mesg);
    for (i=0; i < (m+2)*(n+2); i++){
       int rowIndex = i / (n+2);
       int colIndex = i % (n+2);
       if ((colIndex>0) && (colIndex<m+1))
          if ((rowIndex > 0) && (rowIndex < n+1))
            printf("%6.3f ", E[i]);
       if (colIndex == m+1)
	    printf("\n");
    }
}
