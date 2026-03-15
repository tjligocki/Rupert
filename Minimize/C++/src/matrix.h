#ifndef MATRIX_H
#define MATRIX_H

class MATRIX
{
  public:
    MATRIX(int a_n, int a_m);

    MATRIX(const MATRIX& a_mat) = delete;
    MATRIX& operator=(const MATRIX& a_mat) = delete;

    ~MATRIX();

    void zero();

    void identity();

    void rotate(double a, int i1, int i2);

    double det();

    void print();

    int m_n,m_m;
    double **m_mat;
};

#endif
