/* 
 * Driver for a cardiac elecrophysioly simulatin that uses the
 * Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD, 1/11/2016
 */

#include <iostream>
#include <assert.h>
#include "apf.h"
#include "Plotting.h"
using namespace std;

// External functions
double *alloc1D(int m,int n);
void cmdLine(int argc, char *argv[]);
void init (double *E,double *E_prev,double *R,int m,int n);
void ReportStart(double dt);
void ReportEnd(double l2norm, double mx, double t0);
double getTime();
void solve(double **_E, double **_E_prev, double *R, double alpha, double dt, Plotter *plotter, double &L2, double &Linf);
void Stop();

control_block cb;
// Main program
int main(int argc, char** argv)
{
 /*
  *  Solution arrays
  *   E is the "Excitation" variable, a voltage
  *   R is the "Recovery" variable
  *   E_prev is the Excitation variable for the previous timestep,
  *      and is used in time integration
  */
 double *E, *R, *E_prev;

 // Default values for the command line arguments
 cb.m = cb.n=255;
 cb.stats_freq = 0;
 cb.plot_freq = 0;
 cb.nt = 1;
 cb.TID = 0;
 cb.niters=100;
 cb.debug = false;

// Parse command line arguments
 cmdLine( argc, argv);
// The algorithm fails when n is too small
 if (cb.n <= 25){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    Stop();
 }
 cb.m = cb.n;
 int nprocs=1, myrank=0;
 if (cb.debug)
     cb.stats_freq = MAX(cb.stats_freq,1);

 // Allocate contiguous memory for solution arrays
 // The computational box is defined on [1:m+1,1:n+1]
 // We pad the arrays in order to facilitate differencing on the 
 // boundaries of the computational box
 E = alloc1D(cb.m+2,cb.n+2);
 E_prev = alloc1D(cb.m+2,cb.n+2);
 R = alloc1D(cb.m+2,cb.n+2);

 init(E,E_prev,R,cb.m,cb.n);

 //
 // Initizialize various parmaters used by the numerical scheme
 //

 // Compute the timestep dt
 // We should always use double precision values for the folowing variables:
 //    rp, dte, dtr, ddt
 //
 // This ensures that the computation of dte and especially dt
 // will not lose precision (i.e. if computed as single precision values)

 double dx = 1.0/(cb.n-1);
 double rp= kk*(b+1)*(b+1)/4;
 double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
 double dtr=1/(epsilon+((M1/M2)*rp));
 double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
 double dt = (double) ddt;
 double alpha = d*dt/(dx*dx);

 // End Initialization

 // Report various information
 ReportStart(dt);

 Plotter *plotter = NULL;
 if (cb.plot_freq){
     plotter = new Plotter();
     assert(plotter);
 }

 // Start the timer
 double t0 = -getTime();
 double L2, Linf;
 solve(&E, &E_prev, R, alpha, dt, plotter,L2,Linf);

 t0 += getTime();

 // Report various information
 ReportEnd(L2,Linf,t0);

 if (cb.plot_freq){
    cout << "\n\nEnter any input to close the program and the plot...";
    int resp;
    cin >> resp;
  }

 free (E);
 free (E_prev);
 free (R);
 if (cb.plot_freq)
     delete plotter;
}
