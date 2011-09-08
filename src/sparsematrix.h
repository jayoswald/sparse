#pragma once
class Vector;
class SymmetricSparseMatrix;

// Minimal sparse matrix class with CRS storage.
class SparseMatrix {
friend Vector *sym_spmv(const SparseMatrix&, const Vector&, int);
friend SymmetricSparseMatrix* make_sparse(int, int);
public:
    SparseMatrix(int rows, int nnz) 
    : _ia(new int[rows+1]), 
      _ja(new int[nnz]),
      _data(new double[nnz]),
      _nnz(nnz), _rows(rows)
    {}

    ~SparseMatrix() {
        delete [] _ia;
        delete [] _ja;
        delete [] _data;
    }

protected:
    double *_data, *_cholesky;
    int    *_ia, *_ja, _nnz, _rows;
};


class SymmetricSparseMatrix : public SparseMatrix {
friend void test_chol();
public:
    SymmetricSparseMatrix(int rows, int nnz) 
    : _incmplt_cholesky(0),
      SparseMatrix(rows, nnz) 
    {}
    
    ~SymmetricSparseMatrix() { 
        delete [] _incmplt_cholesky;
    }
    // Builds the incomplete Cholesky factorization.
    void init_chol_preconditioner();

protected:
    // Incomplete Cholesky factorization (for preconditioning).
    double *_incmplt_cholesky;

};

// Computes y = A*x for symmetric A. 
// Returns a pointer to y that should be deleted elsewhere.
Vector* sym_spmv(const SparseMatrix &A, const Vector &x, int nthreads);

// Makes a symmetric sparse matrix with nc nonzeros per column (diagonal).
SymmetricSparseMatrix* make_sparse(int n, int bw);

void test_chol();

