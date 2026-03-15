#include <cstdio>
#include <cmath>

#include "matrix.h"

MATRIX::MATRIX(int a_n, int a_m)
{
  m_mat = new double* [a_n];

  for (int i = 0; i < a_n; i++)
  {
    m_mat[i] = new double [a_m];
  }

  m_n = a_n;
  m_m = a_m;
}

MATRIX::~MATRIX()
{
  for (int i = 0; i < m_n; i++)
  {
    delete[] m_mat[i];
  }

  delete[] m_mat;
}

void MATRIX::zero()
{
  for (int i = 0; i < m_n; i++)
  {
    for (int j = 0; j < m_m; j++)
    {
      m_mat[i][j] = 0.0;
    }
  }
}

void MATRIX::identity()
{
  for (int i = 0; i < m_n; i++)
  {
    for (int j = 0; j < m_m; j++)
    {
      if (i != j)
      {
        m_mat[i][j] = 0.0;
      } else {
        m_mat[i][j] = 1.0;
      }
    }
  }
}

void MATRIX::rotate(double a_a, int a_i1, int a_i2)
{
  double sa,ca;

  double **scan;

  scan = m_mat;

  sa = sin(a_a);
  ca = cos(a_a);

  for (int i = 0; i < m_n; i++)
  {
    double e1,e2;

    e1 = ca * scan[a_i1][i] - sa * scan[a_i2][i];
    e2 = sa * scan[a_i1][i] + ca * scan[a_i2][i];

    scan[a_i1][i] = e1;
    scan[a_i2][i] = e2;
  }
}

double MATRIX::det()
{
  if (m_n != m_m)
  {
    fprintf(stderr,"MATRIX::det: matrix dimensions (%d x %d) incompatible\n",m_n,m_m);
    exit(1);
  }

  double **scan = m_mat;
  double sum = 0.0;

  if (m_n == 2)
  {
    sum = scan[0][0] * scan[1][1] - scan[1][0] * scan[0][1];
  } else {
    double parity;
    MATRIX reduced(m_n-1,m_m-1);
    double **scan_r;

    parity = 1.0;
    scan_r = reduced.m_mat;

    for (int i = 0; i < m_n; i++)
    {
      for (int j = 1; j < m_n; j++)
      {
        for (int k = 0; k < m_m; k++)
        {
          if (k < i)
          {
            scan_r[j-1][k] = scan[j][k];
          } else if (k > i) {
            scan_r[j-1][k-1] = scan[j][k];
          }
        }
      }

      sum += parity * scan[0][i] * reduced.det();

      parity *= -1.0;
    }
  }

  return sum;
}

void MATRIX::print()
{
  double **scan = m_mat;

  for (int i = 0; i < m_n; i++)
  {
    for (int j = 0; j < m_m; j++)
    {
      printf("%9.6f ",scan[i][j]);
    }

    printf("\n");
  }
}
