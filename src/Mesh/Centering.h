// *****************************************************************************
/*!
  \file      src/Mesh/Centering.h
  \copyright 2016-2018, Triad National Security, LLC.
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
