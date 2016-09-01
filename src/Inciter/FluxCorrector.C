// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.C
  \author    J. Bakosi
  \date      Thu 01 Sep 2016 08:22:19 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. There is a
    FluxCorrector Charm++ array element bound to each Carrier array element..
    Each FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************

#include "Vector.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::FluxCorrector;

FluxCorrector::FluxCorrector()
// *****************************************************************************
//  Constructor
//! \author J. Bakosi
// *****************************************************************************
{
}

void
FluxCorrector::fct( const std::array< std::vector< tk::real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const tk::MeshNodes& u,
                    const tk::MeshNodes& uh )
// *****************************************************************************
//  Start FCT
//! \author J. Bakosi
// *****************************************************************************
{
  aec( coord, inpoel, u, uh );
}

void
FluxCorrector::aec( const std::array< std::vector< tk::real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const tk::MeshNodes& Un,
                    const tk::MeshNodes& Uh )
// *****************************************************************************
//  Compute antidiffusive element contributions (AEC)
//! \author J. Bakosi
// *****************************************************************************
{
  auto dUh = Uh - Un;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

  std::vector< tk::real > AEC( inpoel.size(), 0.0 );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );

    // construct tetrahedron element-level matrices

    // consistent - lumped mass
    std::array< std::array< tk::real, 4 >, 4 > m;       // nnode*nnode [4][4]
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = J/40.0;     // diagonal
    m[0][1] = m[0][2] = m[0][3] =                       // off-diagonal
    m[1][0] = m[1][2] = m[1][3] =
    m[2][0] = m[2][1] = m[2][3] =
    m[3][0] = m[3][1] = m[3][2] = -J/120.0;

    // access solution at element nodes at time n
    std::vector< std::array< tk::real, 4 > > un( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, 0, N );
    // access high-order solution increment at element nodes
    std::vector< std::array< tk::real, 4 > > duh( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) duh[c] = dUh.extract( c, 0, N );

    // compute anti-diffusive element contributions
    for (ncomp_t c=0; c<ncomp; ++c)
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k)
          AEC[e*4+j] = m[j][k] * (ctau*un[c][k] + duh[c][k]);
  }
}
