#pragma once 

class SparseMatrix;

// Very simple vector class.
class Vector {
friend Vector *sym_spmv(const SparseMatrix&, const Vector&, int);
public:
    Vector(int len) : _data(new double[len]), _len(len) {}
    ~Vector() { delete _data; }

    double &operator()(int i) { return _data[i]; }
    double  operator()(int i) const { return _data[i]; }
    int size() const { return _len; }
private:
    double *_data;
    int     _len;
};

