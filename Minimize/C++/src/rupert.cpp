#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rupert.h"

RUPERT::RUPERT(const ROTATE &rot)
{
  double *scan;

  int i,j;

  m_m = rot.m_m;
  m_n = rot.m_n;

  m_mat = new MATRIX(m_n,m_n);
  m_mat->identity();

  scan = rot.m_rot + rot.m_total;
  for (int i = m_m-1; i >= 0; i--)
  {
    for (int j = m_n-1; j > i; j--)
    {
      scan--;
      m_mat->rotate(*scan,i,j);
    }
  }
}

double RUPERT::value()
{
  double **scan = m_mat->m_mat;

  double r = 0.0;

  for (int i = 0; i < m_n; i++)
  {
    double sum = 0.0;

    for (int j = 0; j < m_m; j++)
    {
      sum += fabs(scan[i][j]);
    }

    if (sum > r)
    {
      r = sum;
    }
  }

  return r;
}

RUPERT::~RUPERT()
{
  delete m_mat;
}
