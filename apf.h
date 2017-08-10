#ifndef apf_h
#define apf_h
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
 // Various constants - these definitions shouldn't change
 const double a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;

#define MAX(a,b) ((a) > (b) ) ? (a) : (b)
#define MIN(a,b) ((a) < (b) ) ? (a) : (b)

typedef struct _control_block {
    int m,n;
    int stats_freq;
    int plot_freq;
    int nt;
    int TID;
    int niters;
    bool debug;
    bool wait;
} control_block;

extern control_block cb;

#else
#endif
