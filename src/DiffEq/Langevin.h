// *****************************************************************************
/*!
  \file      src/DiffEq/Langevin.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Functionality implementing Langevin models for the velocity
  \details   Functionality implementing Langevin models for the velocity.
*/
// *****************************************************************************
#ifndef Langevin_h
#define Langevin_h

#include <array>

#include "Types.h"
#include "Walker/Options/Depvar.h"

namespace walker {

//! Calculate the 2nd order tensor Gij based on the simplified Langevin model
std::array< tk::real, 9 >
slm( tk::real hts, tk::real C0 );

//! Calculate the 2nd order tensor Gij based on the generalized Langevin model
std::array< tk::real, 9 >
glm( tk::real hts,
     tk::real C0,
     const std::array< tk::real, 6 >& rs,
     const std::array< tk::real, 9 >& dU );

//! Compute the Reynolds stress tensor
std::array< tk::real, 6 >
reynoldsStress( char depvar,
                ctr::DepvarType solve,
                const std::map< tk::ctr::Product, tk::real >& moments );

//! Compute the turbulent kinetic energy
tk::real
tke( char depvar,
     ctr::DepvarType solve,
     const std::map< tk::ctr::Product, tk::real >& moments );

} // walker::

#endif // Langevin_h
