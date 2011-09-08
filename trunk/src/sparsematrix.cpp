#include "sparsematrix.h"
#include <iostream>
#include <thread>
#include <vector>
#include "vector.h"
#include <cmath>

using std::thread;
using std::mutex;
using std::cout;

// Computes y = A*x for symmetric A. 
// Returns a pointer to y that should be deleted elsewhere.
Vector* sym_spmv(const SparseMatrix &A, const Vector &x, int nthreads) {

    std::vector<Vector*> rhs(nthreads);   
    std::vector<std::thread> pool;
    std::vector<mutex*> m(nthreads);

    for (int t=0; t<nthreads; ++t) {
        m[t] = new mutex();
        m[t]->lock();
    }

    for (int t=0; t<nthreads; ++t) {
        pool.push_back(thread(
        // Lambda expression, t captured by value
        // everything else is by reference.
        [&, t]()  {
            int n = x.size();
            int nthreads = rhs.size();
            
            // Create a container and zero it.
            Vector *r = new Vector(n);
            for (int i=0; i<n; ++i) (*r)(i) = 0.0;
       
            // Determine what portion of the job this thread does. 
            int chunk = n / nthreads; 
            int begin = t*chunk;
            int end   = begin + chunk;
            if (t == nthreads-1) end = n; 
            for (int i=begin; i<end; ++i) { 
                for (int ij=A._ia[i]; ij<A._ia[i+1]; ++ij) {
                    int j = A._ja[ij];
                    (*r)(i) += A._data[ij] * x(j);
                    if (i == j) continue;
                    (*r)(j) += A._data[ij] * x(i);
                }
            }
            rhs[t] = r;
            m[t]->unlock(); 
            for (int i=0; i<m.size(); ++i) {
                m[i]->lock();
                m[i]->unlock();
            }
            
            // Copy segment that this thread owns, [beg, end), to r0.
            for (int i=begin; i<end; ++i) {
                for (int td=1; td<nthreads; ++td) {
                    (*rhs[0])(i) += (*rhs[td])(i);
                }
            }
            delete m[t];
        }));
    }
    for (auto i=pool.begin(); i!=pool.end(); ++i) i->join();    
    return rhs[0];
}


// Builds an imcomplete Cholesky factorization.
void SymmetricSparseMatrix::init_chol_preconditioner() {
    if (!_incmplt_cholesky) _incmplt_cholesky = new double[_nnz];
    double *c = _incmplt_cholesky;
    for (int i=0; i<_nnz; ++i)  c[i] = _data[i];
    
    for (int i=0; i<_rows; ++i) {
        double Aii_inv = 0.0;
        // Find diagonal and sqrt it.
        for (int ij=_ia[i]; ij!=_ia[i+1]; ++ij) {
            if (i == _ja[ij]) {
                c[ij] = sqrt(c[ij]);
                Aii_inv = 1.0/c[ij]; 
                break;
            }
        }
        // Divides L by row by sqrt(diagonal.)
        for (int ij=_ia[i]; ij!=_ia[i+1]; ++ij) {
            if (i != _ja[ij]) c[ij]*=Aii_inv;
        }

        for (int ij=_ia[i]; ij!=_ia[i+1]; ++ij) {
            int j = _ja[ij];
            if (j == i) continue;
            for (int ik=_ia[i]; ik!=_ia[i+1]; ik++) {
                if (_ja[ik] == i) continue;
                // Find k in row j.
                for (int jk=_ia[j]; jk!=_ia[j+1]; ++jk) {
                    if (_ja[jk] == _ja[ik]) {
                        c[jk] -= c[ik]*c[ij];
                        break;
                    }
                }
            }
        }
    }
}

// Makes a sparse matrix with nc nonzeros per column (diagonal).
SymmetricSparseMatrix *make_sparse(int n, int bw) {
    int nnz = bw * n; 
    SymmetricSparseMatrix *A = new SymmetricSparseMatrix(n, nnz);
    nnz = 0;
    for (int i=0; i<n; ++i) {
        A->_ia[i] = nnz;
        for (int j=i; j<std::min(n, i+bw); ++j) {

            A->_ja[nnz] = j;
            A->_data[nnz]  = 2.0*double(i==j) - 0.5;
            nnz++;
        }        
    }
    A->_ia[n] = nnz;
    std::cout << "Created a matrix with " << 2*nnz -n << " nonzeros.\n";
    return A;
}

void test_chol() {
    int ia[]   = {0,          4,       7,    9, 10};
    int ja[]   = {0, 1, 2, 3, 1, 2, 3, 2, 3, 3};
    double v[] = {2, 1, 1, 1, 2, 1, 1, 2, 1, 2};

    int rows = sizeof(ia)/sizeof(int)-1;
    int nnz = sizeof(v)/sizeof(double);
    SymmetricSparseMatrix A(rows, nnz);
    for (int i=0; i<rows+1; ++i) A._ia[i] = ia[i];

    for (int i=0; i<nnz; ++i) {
        A._ja[i]   = ja[i];
        A._data[i] = v[i];
    }
    A.init_chol_preconditioner();

    for (int i=0; i<nnz; ++i) {
        cout << A._incmplt_cholesky[i] << "\n";
    }

}

