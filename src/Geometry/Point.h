//******************************************************************************
/*!
  \file      src/Geometry/Point.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:04:50 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Point primitive
  \details   Point primitive
*/
//******************************************************************************
#ifndef Point_h
#define Point_h

#include <QuinoaTypes.h>

namespace quinoa {

//! Point in 3D space
struct Point {
  real x;
  real y;
  real z;

  //! Fill constructor
  explicit Point(const real X, const real Y, const real Z) : x(X), y(Y), z(Z) {}
};

} // namespace quinoa

#endif // Point_h
