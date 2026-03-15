#ifndef ANGLE_H
#define ANGLE_H

inline double rad2deg(double a_rad)
{
  return a_rad / M_PI * 180.0;
}

inline double deg2rad(double a_deg)
{
  return a_deg / 180.0 * M_PI;
}

#endif
