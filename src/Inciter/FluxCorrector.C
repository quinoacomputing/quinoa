// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.C
  \author    J. Bakosi
  \date      Mon 12 Sep 2016 03:55:47 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. There is a
    FluxCorrector Charm++ array element bound to each Carrier array element..
    Each FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************

#include <limits>
#include <algorithm>

#include "Vector.h"
#include "DerivedData.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

using inciter::FluxCorrector;

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
  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  Assert( coord[0].size() == coord[1].size() &&
          coord[0].size() == coord[2].size(),
          "Mesh node coordinate size mismatch" );
  Assert( coord[0].size() == Un.nunk() && Un.nunk() == Uh.nunk(),
          "Mesh node coordinate and unknown array size mismatch" );
  Assert( m_aec.nunk() == inpoel.size() && m_aec.nprop() == ncomp,
          "AEC and mesh connectivity size mismatch" );

  auto dUh = Uh - Un;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

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
          m_aec(e*4+j,c,0) = m[j][k] * (ctau*un[c][k] + duh[c][k]);
  }

  tk::MeshNodes P( Un.nunk(), Un.nprop()*2 );
  P.fill( 0.0 );

  // sum all positive (negative) antidiffusive element contributions to nodes
  // (Lohner: P^{+,-}_i)
  const auto esup = tk::genEsup( inpoel, 4 );
  for (std::size_t p=0; p<Un.nunk(); ++p)
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i) {
       const auto e = esup.first[i];
       // decide which node's contribution we need (A, B, C, or D)
       std::size_t n = 3;
            if (inpoel[e*4+0] == p) n = 0;
       else if (inpoel[e*4+1] == p) n = 1;
       else if (inpoel[e*4+2] == p) n = 2;
       // add up element contributions to node p
       for (ncomp_t c=0; c<ncomp; ++c) {
         P(p,c,0) += std::max( 0.0, m_aec(e*4+n,c,0) );
         P(p,c,1) += std::min( 0.0, m_aec(e*4+n,c,0) );
       }
    }

  return P;
}

tk::MeshNodes
FluxCorrector::lump( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel ) const
// *****************************************************************************
//  Compute lumped mass matrix lhs for low order system
//! \return Lumped mass matrix
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  tk::MeshNodes L( coord[0].size(), ncomp );
  L.fill( 0.0 );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) * 5.0 / 120.0;

    // access pointer to lumped mass left hand side at element nodes
    std::vector< const tk::real* > l( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) l[c] = L.cptr( c, 0 );

    // scatter-add lumped mass element contributions to lhs nodes
    for (ncomp_t c=0; c<ncomp; ++c)
      for (std::size_t j=0; j<4; ++j)
        L.var(l[c],N[j]) += J;
  }

  return L;
}

tk::MeshNodes
FluxCorrector::diff( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const tk::MeshNodes& Un ) const
// *****************************************************************************
//  Compute mass diffusion contribution to the rhs of the low order system
//! \return Mass diffusion contribution to the rhs of the low order system
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

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );

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

    // scatter-add mass diffusion element contributions to rhs nodes
    for (ncomp_t c=0; c<ncomp; ++c)
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k)
           D.var(d[c],N[j]) -= ctau * m[j][k] * un[c][k];
  }

  return D;
}

tk::MeshNodes
FluxCorrector::allowed( const std::vector< std::size_t >& inpoel,
                        const tk::MeshNodes& Un,
                        const tk::MeshNodes& Ul) const
// *****************************************************************************
//  Compute the allowed solution increments and decrements at mesh nodes
//! \author J. Bakosi
// *****************************************************************************
{
  // compute maximum and minimum nodal values of Ul and Un (Lohner: u^*_i)
  auto Smax = tk::max( Ul, Un );
  auto Smin = tk::min( Ul, Un );

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  // compute maximum and minimum nodal values of all elements (Lohner: u^*_el)
  tk::MeshNodes S( inpoel.size()/4, ncomp*2 );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      S(e,c,0) = std::max( Smax(N[3],c,0),
                   std::max( Smax(N[2],c,0),
                     std::max( Smax(N[0],c,0), Smax(N[1],c,0) ) ) );
      S(e,c,1) = std::min( Smin(N[3],c,0),
                   std::min( Smin(N[2],c,0),
                     std::min( Smin(N[0],c,0), Smin(N[1],c,0) ) ) );
    }
  }

  // compute maximum and mimimum unknowns of all elements surrounding each node
  // (Lohner: u^{max,min}_i)
  tk::MeshNodes Q( Un.nunk(), Un.nprop()*2 );
  const auto esup = tk::genEsup( inpoel, 4 );
  for (std::size_t p=0; p<Un.nunk(); ++p) {
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c,0) = -std::numeric_limits< tk::real >::max();
      Q(p,c,1) = std::numeric_limits< tk::real >::max();
    }
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i) {
      const auto e = esup.first[i];
      for (ncomp_t c=0; c<ncomp; ++c) {
        if (S(e,c,0) > Q(p,c,0)) Q(p,c,0) = S(e,c,0);
        if (S(e,c,1) < Q(p,c,1)) Q(p,c,1) = S(e,c,1);
      }
    }
  }

  // compute the maximum and minimum increments and decrements nodal solution
  // values are allowed to achieve (Lohner: Q^{+,-}_i)
  for (std::size_t p=0; p<Ul.nunk(); ++p)
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c,0) -= Ul(p,c,0);
      Q(p,c,1) -= Ul(p,c,0);
    }

  // maximum and minimum increments and decrements nodal solution values are
  // allowed to achieve
  return Q;
}

void
FluxCorrector::limit( const std::vector< std::size_t >& inpoel,
                      const tk::MeshNodes& P,
                      tk::MeshNodes& Q,
                      tk::MeshNodes& U ) const
// *****************************************************************************
// Perform limiting
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( P.nunk() == Q.nunk() && P.nprop() == Q.nprop(), "Size mismatch" );
  Assert( P.nunk() == U.nunk() && P.nprop() == U.nprop()*2, "Size mismatch" );

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  // compute the ratios of positive and negative element contributions that
  // ensure monotonicity (Lohner: R^{+,-})
  for (std::size_t p=0; p<P.nunk(); ++p)
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c,0) = (P(p,c,0) > 0.0 ? std::min(1.0,Q(p,c,0)/P(p,c,0)) : 0.0);
      Q(p,c,1) = (P(p,c,1) < 0.0 ? std::min(1.0,Q(p,c,1)/P(p,c,1)) : 0.0);
    }

  // calculate limit coefficient for all elements (Lohner: C_el)
  tk::MeshNodes C( inpoel.size()/4, ncomp );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      std::array< tk::real, 4 > R;
      for (std::size_t j=0; j<4; ++j)
        R[j] = m_aec(e*4+j,c,0) > 0.0 ? Q(N[j],c,0) : Q(N[j],c,1);
      C(e,c,0) = *std::min_element( begin(R), end(R) );
      Assert( C(e,c,0) > -std::numeric_limits< tk::real >::epsilon() &&
              C(e,c,0) < 1.0+std::numeric_limits< tk::real >::epsilon(),
              "0 <= AEC <= 1.0 failed" );
    }
  }

  // apply the limited antidiffusive element contributions
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};

    // access pointer to solution at element nodes
    std::vector< const tk::real* > u( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) u[c] = U.cptr( c, 0 );

    // scatter-add limited antidiffusive element contributions to nodes
    for (ncomp_t c=0; c<ncomp; ++c) 
      for (std::size_t j=0; j<4; ++j)
        U.var(u[c],N[j]) += C(e,c,0) * m_aec(e*4+j,c,0);
  }
}
