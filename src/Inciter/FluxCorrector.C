// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.C
  \author    J. Bakosi
  \date      Tue 06 Sep 2016 12:27:30 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. There is a
    FluxCorrector Charm++ array element bound to each Carrier array element..
    Each FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************

#include <algorithm>

#include "Vector.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::FluxCorrector;

FluxCorrector::FluxCorrector( std::pair< std::vector< std::size_t >,
                                         std::vector< std::size_t > >&& esup )
  : m_esup( std::move(esup) )
// *****************************************************************************
//  Constructor
//! \author J. Bakosi
// *****************************************************************************
{
}

tk::MeshNodes
FluxCorrector::aec( const std::array< std::vector< tk::real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const tk::MeshNodes& Un,
                    const tk::MeshNodes& Uh )
// *****************************************************************************
//  Compute antidiffusive element contributions (AEC)
//! \details The antidiffusive element contributions (AEC) are computed as the
//!   difference between the high and low order solution, where the high order
//!   solution is consistent mass Taylor-Galerkin and the low order solution is
//!   lumped mass Taylor-Galerkin + diffusion.
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( coord[0].size() == coord[1].size() &&
          coord[0].size() == coord[2].size(),
          "Mesh node coordinate size mismatch" );
  Assert( coord[0].size() == Un.nunk() && Un.nunk() == Uh.nunk(),
          "Mesh node coordinate and unknown array size mismatch" );

  auto dUh = Uh - Un;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

  tk::MeshNodes AEC( inpoel.size(), ncomp );

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

    // lumped - consistent mass
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

    // compute antidiffusive element contributions
    for (ncomp_t c=0; c<ncomp; ++c)
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k)
          AEC(e*4+j,c,0) = m[j][k] * (ctau*un[c][k] + duh[c][k]);
  }

  tk::MeshNodes P( Un.nunk(), Un.nprop()*2 );
  P.fill( 0.0 );

  // sum all positive (negative) antidiffusive element contributions to nodes
  for (std::size_t p=0; p<Un.nunk(); ++p)
    for (auto i=m_esup.second[p]+1; i<=m_esup.second[p+1]; ++i) {
       auto e = m_esup.first[i];
       // decide which node's contribution we need (A, B, C, or D)
       std::size_t n = 3;
            if (inpoel[e*4+0] == p) n = 0;
       else if (inpoel[e*4+1] == p) n = 1;
       else if (inpoel[e*4+2] == p) n = 2;
       // add up element contributions to node p
       for (ncomp_t c=0; c<ncomp; ++c) {
         P(p,c*2+0,0) += std::max( 0.0, AEC(e*4+n,c,0) );
         P(p,c*2+1,0) += std::min( 0.0, AEC(e*4+n,c,0) );
       }
    }

  return P;
}

std::pair< tk::MeshNodes, tk::MeshNodes >
FluxCorrector::low( const std::array< std::vector< tk::real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const tk::MeshNodes& Un )
// *****************************************************************************
//  Compute lhs and rhs for low order solution
//! \details Compute lumped mass for as the left-hand side, and mass diffusion
//!   to be added to the low order solution rhs
//! \return Lumped mass matrix and mass diffusion contribution of the rhs
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

  tk::MeshNodes D( Un.nunk(), Un.nprop() );
  D.fill( 0.0 );
  tk::MeshNodes L( Un.nunk(), Un.nprop() );
  L.fill( 0.0 );

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

    // lumped - consistent mass
    std::array< std::array< tk::real, 4 >, 4 > m;       // nnode*nnode [4][4]
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = J/40.0;     // diagonal
    m[0][1] = m[0][2] = m[0][3] =                       // off-diagonal
    m[1][0] = m[1][2] = m[1][3] =
    m[2][0] = m[2][1] = m[2][3] =
    m[3][0] = m[3][1] = m[3][2] = -J/120.0;

    // access solution at element nodes at time n
    std::vector< std::array< tk::real, 4 > > un( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, 0, N );
    // access pointer to mass diffusion right hand side at element nodes
    std::vector< const tk::real* > d( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) d[c] = D.cptr( c, 0 );
    // access pointer to lumped mass left hand side at element nodes
    std::vector< const tk::real* > l( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) l[c] = L.cptr( c, 0 );

    for (ncomp_t c=0; c<ncomp; ++c) 
      for (std::size_t j=0; j<4; ++j) {
        // scatter-add lumped mass element contributions to lhs nodes
        L.var(l[c],N[j]) += 5.0*J/120.0;
        // scatter-add mass diffusion element contributions to rhs nodes
        for (std::size_t k=0; k<4; ++k)
           D.var(d[c],N[j]) -= ctau * m[j][k] * un[c][k];
      }
  }

  return { L, D };
}
