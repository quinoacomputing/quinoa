// *****************************************************************************
/*!
  \file      src/PDE/Integrate/RhsAccum.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Accumulators for the const-order (_constP) rhs updates
  \details   Overload pairs that let one templated body write either to host
             tk::Fields / std::vector or to device Kokkos::Views. The device
             overloads are atomic, which internal faces require: a face is
             written by both of its adjacent elements, so a face-parallel
             kernel has two writers per location. Boundary faces write to one
             element only and do not need atomicity, but use the same
             overloads for uniformity.

             Kept in a header rather than an anonymous namespace so that both
             Surface.cpp and Boundary.cpp can use one definition.
*/
// *****************************************************************************
#ifndef RhsAccum_h
#define RhsAccum_h

#include <vector>

#include "Fields.hpp"
#include "Types.hpp"
#include "KokkosDevice.hpp"

namespace tk {

//! Accumulate into R, host tk::Fields
inline void
rhsAccum( tk::Fields& R, std::size_t e, std::size_t /*nprop*/,
          std::size_t idx, tk::real v )
{
  R(e,idx) += v;
}

//! Accumulate into R, Kokkos::View (device, atomic)
KOKKOS_INLINE_FUNCTION
void
rhsAccum( const Kokkos::View< tk::real*, memory_space >& R, std::size_t e,
          std::size_t nprop, std::size_t idx, tk::real v )
{
  Kokkos::atomic_add( &R(e*nprop + idx), v );
}

//! Accumulate into riemannDeriv, host vector of vectors
inline void
rdAccum( std::vector< std::vector< tk::real > >& rd, std::size_t row,
         std::size_t col, std::size_t /*ncol*/, tk::real v )
{
  rd[row][col] += v;
}

//! Accumulate into riemannDeriv, Kokkos::View (device, atomic)
KOKKOS_INLINE_FUNCTION
void
rdAccum( const Kokkos::View< tk::real*, memory_space >& rd, std::size_t row,
         std::size_t col, std::size_t ncol, tk::real v )
{
  Kokkos::atomic_add( &rd(row*ncol + col), v );
}

} // tk::

#endif // RhsAccum_h
