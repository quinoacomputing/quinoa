// *****************************************************************************
/*!
  \file      src/Inciter/ALE/MeshMotion.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh motion functionality for ALE
  \details   Functionality for arbitrary Lagrangian-Eulerian (ALE) mesh motion.
*/
// *****************************************************************************
#ifndef MeshMotion_h
#define MeshMotion_h

#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "Control/Inciter/Options/MeshVelocity.hpp"

namespace inciter {

//! Prescribe mesh motion as a sine wave for testing purposes
void sine( const tk::UnsMesh::Coords& coord, tk::Fields& w );

//! Assign mesh velocity based on user config
void
meshvel( ctr::MeshVelocityType m,
         const tk::UnsMesh::Coords& coord,
         const tk::UnsMesh::Coords& v,
         tk::Fields& w );

//! Scale the mesh velocity with a function of the fluid vorticity for ALE
void
vortscale( const std::vector< tk::real >& vort,
           tk::real vmult,
           tk::real maxv,
           tk::Fields& w );

} // inciter::

#endif // MeshMotion_h
