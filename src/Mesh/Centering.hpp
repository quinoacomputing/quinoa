// *****************************************************************************
/*!
  \file      src/Mesh/Centering.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh solution location (centering)
  \details   Mesh solution location (centering) enum type.
*/
// *****************************************************************************
#ifndef Centering_h
#define Centering_h

namespace tk {

//! Mesh/scheme centering types
//! \see Control/Inciter/Options/Scheme.h
enum class Centering : uint8_t { NODE
                               , ELEM };

} // tk::

#endif // Centering_h
