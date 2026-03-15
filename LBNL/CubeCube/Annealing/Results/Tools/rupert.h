#ifndef _rupert_h_
#define _rupert_h_

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>
using std::vector;

#include <mpfr.h>

#define NUMBER_PRECISION  300

class Number {
  public:
    Number();

    Number(const Number& in);

    ~Number();

    const Number& operator=(const Number& in);

    inline mpfr_ptr Get() {
      return number;
    };

  private:
    mpfr_t number;
};

class Rupert {
  public:
    Rupert(int m0, int n0);

    Rupert(const Rupert& in);

    ~Rupert();

    const Rupert& operator=(const Rupert& in);

    void Read();

    void Write();

    void DotRows(Number& dot, int i1, int i2);

    void Compute();

    void OrthoNormal();

    void Gradient(Number& grad, int i1, int i2, double dt);

    void Transform(Rupert& trans, int i1, int i2, const Number& angle);

    vector<vector<Number> >& GetMatrix();

    vector<Number>& GetSum();

    Number& GetMaxSum();

    Number& GetSide();

  private:
    vector<vector<Number> > matrix;
    vector<Number>          sum;
    Number                  maxSum;
    Number                  side;
    int                     m;
    int                     n;
};

#endif
