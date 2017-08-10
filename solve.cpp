/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <thread>
#include <math.h>
#include "time.h"
#include "apf.h"
#include <mutex>
#include "Plotting.h"

#ifdef SSE_VEC
#include <emmintrin.h>
#endif 
using namespace std;

//barrier class definition

class Barrier{
   int count, _NT;
   mutex arrival, departure;

  public:
   Barrier(int NT){
      count = 0;
      _NT = NT;
      arrival.unlock();
      departure.lock();
   }

   void bsync(){
      arrival.lock();
      if(++count < _NT){
        arrival.unlock();}
      else 
        departure.unlock();
      departure.lock();
      if(--count > 0){
        departure.unlock();}
      else
        arrival.unlock();
   }
};

void repNorms(double l2norm, double mx, double dt, int m,int n, int niter, int stats_freq);
void stats(double *E, int m, int n, double *_mx, double *sumSq);

#ifdef SSE_VEC
// If you intend to vectorize using SSE instructions, you must
// disable the compiler's auto-vectorizer
__attribute__((optimize("no-tree-vectorize")))
#endif 


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
double L2Norm(double sumSq){
    double l2norm = sumSq /  (double) ((cb.m)*(cb.n));
    l2norm = sqrt(l2norm);
    return l2norm;
}


//thread function 
void threadfunc(int start, int end, Barrier* bar, double** _E, double ** _E_prev, 
    double *R, double alpha, double dt, Plotter *plotter, double L2, 
    double Linf, int m, int n, double mx, double sumSq, double * E, 
    double * E_prev,int id)
{
 double t = 0.0;
 double *R_tmp = R;
 double *E_tmp = *_E;
 double *E_prev_tmp = *_E_prev;
 int niter;

 // We continue to sweep over the mesh until the simulation has reached
 // the desired number of iterations
  for (niter = 0; niter < cb.niters; niter++){
      if  (cb.debug && (niter==0)){
	      stats(E_prev,m,n,&mx,&sumSq);
         double l2norm = L2Norm(sumSq);
	      repNorms(l2norm,mx,dt,m,n,-1, cb.stats_freq);
          if (cb.plot_freq){
              bar->bsync();
	         plotter->updatePlot(E,  -1, m+1, n+1);
          }
      }

   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    * These are physical boundary conditions, and are not to be confused
    * with ghost cells that we would use in an MPI implementation
    *
    * The reason why we copy boundary conditions is to avoid
    * computing single sided differences at the boundaries
    * which increase the running time of solve()
    *
    */
    
    // 4 FOR LOOPS set up the padding needed for the boundary conditions
    int i,j;
    
      if(id == 0){
    // Fills in the TOP Ghost Cells
    for (i = 0; i < (n+2); i++) {
        E_prev[i] = E_prev[i + (n+2)*2];
    }

    // Fills in the RIGHT Ghost Cells
      for (i = (n+1); i < (m+2)*(n+2); i+=(n+2)) {
        E_prev[i] = E_prev[i-2];
    }

    // Fills in the LEFT Ghost Cells
      for (i = 0; i < (m+2)*(n+2); i+=(n+2)) {
        E_prev[i] = E_prev[i+2];
    }	

    // Fills in the BOTTOM Ghost Cells
    for (i = ((m+2)*(n+2)-(n+2)); i < (m+2)*(n+2); i++) {
        E_prev[i] = E_prev[i - (n+2)*2];
    }
      }
    // Solve for the excitation, a PDE
      for(j = start; j < end; j+=(n+2)) {
        E_tmp = E + j;
        E_prev_tmp = E_prev + j;
        R_tmp = R + j;
        for(int i = 0;i<n;i++){
         E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
        }
        for(int i = 0; i < n; i++){
         E_tmp[i] += -dt*(kk*E_tmp[i]*(E_tmp[i]-a)*(E_tmp[i]-1)+E_tmp[i]*R_tmp[i]);
        }
        for(int i = 0; i < n; i++){
         R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_tmp[i]+M2))*(-R_tmp[i]-kk*E_tmp[i]*(E_tmp[i]-b-1));
        }
      }
      bar->bsync();
   if (cb.stats_freq){
     if ( !(niter % cb.stats_freq)){
        stats(E,m,n,&mx,&sumSq);
        double l2norm = L2Norm(sumSq);
        repNorms(l2norm,mx,dt,m,n,niter, cb.stats_freq);
    }
   }
      /*
      if (cb.plot_freq){
          if (!(niter % cb.plot_freq)){
              bar->bsync();
              mutex mu;
              mu.lock();
              plotter->updatePlot(E,niter, m, n);
              mu.unlock();
          }
      }
       */
   // Swap current and previous meshes
    mutex curr;
    curr.lock();
    double *tmp = E; E = E_prev; E_prev = tmp;
    curr.unlock();
 } //end of 'niter' loop at the beginning
}
//end of thread function


void solve(double **_E, double **_E_prev, double *R, double alpha, double dt, Plotter *plotter, double &L2, double &Linf){

 // Simulated time is different from the integer timestep number
 int m = cb.m, n=cb.n;
 double mx, sumSq;
 double *E = *_E, *E_prev = *_E_prev;
 int chunk = cb.m/cb.nt; 
 int innerBlockRowStartIndex = (n+2)+1;
 int innerBlockRowEndIndex = 0;
 Barrier *bar = new Barrier(cb.nt);
 thread *thrd = new thread[cb.nt-1];
 
 if(cb.nt == 1){
   threadfunc(n+2 + 1, n+2 + 1 + (chunk) *(n+2), bar, _E, _E_prev, R, alpha, dt, 
      plotter, L2, Linf, m, n, mx, sumSq, E, E_prev,0);
 }
 else{
   for(int i = 0; i < cb.nt; i++){
     if(i == cb.nt - 1)
      innerBlockRowEndIndex = (m+2)*(n+2) - n - 1 - n - 1;
     else 
      innerBlockRowEndIndex = innerBlockRowStartIndex + chunk*(n+2);
      //cout << "start is " << innerBlockRowStartIndex <<" end is " << innerBlockRowEndIndex <<endl;
     thrd[i] = thread(threadfunc, innerBlockRowStartIndex, innerBlockRowEndIndex, 
         bar, _E, _E_prev, R, alpha, dt, plotter, L2, Linf, m, n, mx, sumSq,
         E, E_prev,i);
     //cout<<"Thread " << i <<" created!" << endl;
     innerBlockRowStartIndex += chunk*(n+2);
   }
   for(int i = 0; i < cb.nt; i++){
      thrd[i].join();
    // cout<<"Thread " << i <<" joined!" << endl;
   }
 }
  // return the L2 and infinity norms via in-out parameters
  stats(E_prev,m,n,&Linf,&sumSq);
  L2 = L2Norm(sumSq);

  // Swap pointers so we can re-use the arrays
  *_E = E;
  *_E_prev = E_prev;
}





/*original solve function
 //#pragma omp parallel for schedule(dynamic, chunk)
 for(i = 0; i < n; i++) {
 E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
 }
 */


/*original solve
 
 for(i = 0; i < n; i++) {
 E_tmp[i] += -dt*(kk*E_tmp[i]*(E_tmp[i]-a)*(E_tmp[i]-1)+E_tmp[i]*R_tmp[i]);
 R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_tmp[i]+M2))*(-R_tmp[i]-kk*E_tmp[i]*(E_tmp[i]-b-1));
 }
 */






