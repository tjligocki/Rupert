#include <cstdio>
#include <cmath>

#include "rotate.h"

ROTATE::ROTATE(int a_m, int a_n)
{
  m_m = a_m;
  m_n = a_n;

  m_total = (m_m * (2*m_n - m_m - 1)) / 2;

  m_rot = new double [m_total];
}

double ROTATE::norm()
{
  double norm2 = 0.0;

  for (int i = 0; i < m_total; i++)
  {
    norm2 += m_rot[i] * m_rot[i];
  }

  return sqrt(norm2);
}

void ROTATE::copy(const ROTATE &a_rot)
{
  if (m_m != a_rot.m_m || m_n != a_rot.m_n)
  {
    fprintf(stderr,"ROTATE::copy - rotation dimension mismatch, (%d x %d) vs. (%d x %d)\n",m_n,m_m,a_rot.m_m,a_rot.m_n);
    exit(1);
  }

  for (int i = 0; i < m_total; i++)
  {
    m_rot[i] = a_rot.m_rot[i];
  }
}

ROTATE::~ROTATE()
{
  delete[] m_rot;
}
