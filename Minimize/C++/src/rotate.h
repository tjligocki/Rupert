#ifndef ROTATE_H
#define ROTATE_H

class ROTATE
{
  public:
    ROTATE(int a_m, int a_n);

    ROTATE(const ROTATE& a_mat) = delete;
    ROTATE& operator=(const ROTATE& a_mat) = delete;

    ~ROTATE();

    double norm();
    void copy(const ROTATE &a_rot);

    int m_m,m_n,m_total;
    double *m_rot;
};

#endif
