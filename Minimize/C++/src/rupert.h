#ifndef RUPERT_H
#define RUPERT_H

#include "rotate.h"
#include "matrix.h"

class RUPERT
{
  public:
    RUPERT(const ROTATE &a_rot);

    RUPERT(const RUPERT& a_rup) = delete;
    RUPERT& operator=(const RUPERT& a_rup) = delete;

    ~RUPERT();

    double value();

    int m_m,m_n;
    MATRIX *m_mat;
};

#endif
