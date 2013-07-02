//******************************************************************************
/*!
  \file      src/Geometry/Point.h
  \author    J. Bakosi
  \date      Tue Jul  2 13:06:21 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Point primitive
  \details   Point primitive
*/
//******************************************************************************
#ifndef Point_h
#define Point_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! Point in 3D space
struct Point {
  real x;
  real y;
  real z;

  //! Fill constructor
  explicit Point(const real X, const real Y, const real Z) : x(X), y(Y), z(Z) {}
};

} // namespace Quinoa

#endif // Point_h
