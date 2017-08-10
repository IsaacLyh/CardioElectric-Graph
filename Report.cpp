// 
// Performs various reporting functions
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "apf.h"
using namespace std;


// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem

// Report statistics periodically
void repNorms(double l2norm, double mx,  double dt, int m,int n, int niter, int stats_freq){

    cout <<      setw(6);
    cout.setf(ios::fixed);
    cout << "iteration = " << niter << ", ";
    cout.unsetf(ios::fixed);
    cout.setf(ios::scientific);
    cout.precision(6);
    cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
}
void printTOD(string mesg)
{
    time_t tim = time(NULL);
    string s = ctime(&tim);
    cout << endl;
    if (mesg.length() ==  0) {
        cout << "Time of day: " << s.substr(0,s.length()-1) << endl;
    }
    else {
        cout << "[" << mesg << "] " ;
        cout << s.substr(0,s.length()-1) << endl;
    }
    cout << endl;
}


// Computes the gigaflops rate

double gflops(int n, int niter, double time){
    int n2 = n*n;
    int64_t updates = (int64_t) n2 * (int64_t) niter;
    int64_t flops = 28 * updates;
    double flop_rate = (double) flops / time;
    return ( flop_rate/1.0e9);
}


void ReportEnd(double l2norm, double mx, double t0){
    printTOD("Simulation ends");    

    double gf = gflops(cb.n, cb.niters, t0);
    cout << "End at";
    cout <<          setw(6);
    cout.setf(ios::fixed);
    cout << " iteration " << cb.niters << endl;
    cout.unsetf(ios::fixed);
    cout.setf(ios::scientific);
    cout.precision(5);
    cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
    cout.unsetf(ios::scientific);
    cout.unsetf(ios::fixed);
    cout.precision(6);
    cout << "Running Time: " << t0 << " sec.";
    cout.precision(3);
    cout << " [" << gf << " GFlop/sec]" << endl << endl;

    cout << "   M x N    NT    #iter    T_p     Gflops        Linf, L2" << endl;
    cout << "@ " << cb.m << " " << cb.n << " ";
    cout.precision(3);
    cout << "   " << cb.nt << " ";
    cout.precision(6);
    cout <<  "     " << cb.niters << " ";
    cout.precision(4);
    cout << " " << t0 << "    "  << gf << "  ";

    cout.unsetf(ios::fixed);
    cout.setf(ios::scientific);
    cout.precision(5);
    cout << "  " << mx << " " << l2norm << endl;
    cout << "  -----" << endl;

}

void ReportStart(double dt){
    printTOD("Simulation begins");
    cout << "dt= " << dt << ", ";
    cout << "# iters = " << cb.niters << endl;
    cout << "m x n = " << cb.m << " x " << cb.n << endl;
    cout << "# threads: " << cb.nt << endl;

#ifdef SSE_VEC
    cout << "Using SSE Intrinsics\n";
#endif
}

void Stop(){
   cout.flush();
   cerr.flush();
   exit(-1);
}
