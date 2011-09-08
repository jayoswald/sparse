#include <iostream>
#include "vector.h"
#include "sparsematrix.h"
#include "timer.h"
using namespace std;

// Main program.
int main(int narg, char **argv) {

    // Default values.
    int num_rows   = 1000000;
    int band_width = 30;
    int iterations = 16;
    int num_thread = 1;
    
    if (narg > 1) num_rows   = atoi(argv[1]);
    if (narg > 2) band_width = atoi(argv[2]);
    if (narg > 3) iterations = atoi(argv[3]);
    if (narg > 4) num_thread = atoi(argv[4]);

    SymmetricSparseMatrix *A = make_sparse(num_rows, band_width);
    Vector x(num_rows);
    for (int i=0; i<num_rows; ++i) x(i) = 1.0;

/*
    Timer timer1; 
    A->init_chol_preconditioner();
    cout << "Chol time is: " << timer1.elapsed() << "\n";
  */  
    Timer timer; 
    for (int i=0; i<iterations; ++i) {
        Vector *b = sym_spmv(*A, x, num_thread);
        delete b;
    }
    double dt = timer.elapsed();
    dt = dt / double(iterations);
    cout << "mv time is: " << dt << "\n";

    int flop = num_rows * (2*band_width - 1);
    cout << "MFLOPS = " << double(flop) / dt *1e-6 << "\n";
    return 0;
}

