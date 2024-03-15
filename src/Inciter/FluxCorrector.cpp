// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. Each
    FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************

#include <limits>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "Macro.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "FluxCorrector.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"

using inciter::FluxCorrector;

void
FluxCorrector::aec(
  const std::array< std::vector< tk::real >, 3 >& coord,
  const std::vector< std::size_t >& inpoel,
  const std::vector< tk::real >& vol,
  const std::unordered_map< std::size_t,
    std::vector< std::pair< bool, tk::real > > >& bcdir,
  const std::unordered_map< int,
    std::unordered_set< std::size_t > >& symbcnodemap,
  const std::unordered_map< int,
    std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
  const tk::Fields& Un,
  tk::Fields& P )
// *****************************************************************************
//  Compute antidiffusive element contributions (AEC)
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] vol Volume associated to mesh nodes
//! \param[in] bcdir Vector of pairs of bool and boundary condition value
//!   associated to local mesh node IDs at which to set Dirichlet boundary
//!   conditions.
//! \param[in] symbcnodemap Unique set of node ids at which to set symmetry BCs
//!   associated to side set ids
//! \param[in] bnorm Face normals in boundary points: key global node id,
//!   value: unit normal, outer key: side set id
//! \param[in] Un Solution at the previous time step
//! \param[in,out] P The sums of positive (negative) AECs to nodes
//! \details The antidiffusive element contributions (AEC) are defined as the
//!   difference between the high and low order solution, where the high order
//!   solution is obtained from consistent mass Taylor-Galerkin discretization
//!   and the low order solution is lumped mass Taylor-Galerkin + diffusion.
//!   Note that AEC is not directly computed as dUh - dUl (although that could
//!   also be done), but as AEC = M_L^{-1} (M_Le - M_ce) (ctau * Un + dUh),
//!   where
//!    * M_ce is the element's consistent mass matrix,
//!    * M_Le is the element's lumped mass matrix,
//!    * ctau is the mass diffusion coefficient on the rhs of the low order
//!      solution, see also FluxCorrector::diff(),
//!    * Un is the solution at the previous time step
//!    * dUh is the increment of the high order solution, and
//!    * M_L^{-1} is the inverse of the assembled lumped mass matrix, i.e., the
//!      volume associated to a mesh node by summing the quarter of the element
//!      volumes surrounding the node. Note that this is the correct node volume
//!      taking into account that some nodes are on chare boundaries.
//! \note Since we use the lumped-mass for the high-order solution, dUh
//!   does not contribute to AEC, as computed above.
//! \see Löhner, R., Morgan, K., Peraire, J. and Vahdati, M. (1987), Finite
//!   element flux-corrected transport (FEM–FCT) for the Euler and Navier–Stokes
//!   equations. Int. J. Numer. Meth. Fluids, 7: 1093–1109.
//!   doi:10.1002/fld.1650071007
// *****************************************************************************
{
  auto ncomp = g_newinputdeck.get< newtag::ncomp >();
  auto ctau = g_newinputdeck.get< newtag::ctau >();

  Assert( vol.size() == coord[0].size(), "Nodal volume vector size mismatch" );
  Assert( m_aec.nunk() == inpoel.size() && m_aec.nprop() == ncomp,
          "AEC and mesh connectivity size mismatch" );
  Assert( Un.nunk() == P.nunk() && Un.nprop() == P.nprop()/2, "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  m_aec.fill( 0.0 );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};

    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );
    Assert( J > 0, "Element Jacobian non-positive" );

    // lumped - consistent mass
    std::array< std::array< tk::real, 4 >, 4 > m;       // nnode*nnode [4][4]
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 3.0*J/120.0;// diagonal
    m[0][1] = m[0][2] = m[0][3] =                       // off-diagonal
    m[1][0] = m[1][2] = m[1][3] =
    m[2][0] = m[2][1] = m[2][3] =
    m[3][0] = m[3][1] = m[3][2] = -J/120.0;

    // access solution at element nodes at time n
    std::vector< std::array< tk::real, 4 > > un( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, N );

    // Compute antidiffusive element contributions (AEC). The high order system
    // is M_c * dUh = r, where M_c is the consistent mass matrix and r is the
    // high order right hand side. The low order system is constructed from the
    // high order one by lumping the consistent mass matrix and adding mass
    // diffusion: M_L * dUl = r + c_tau * (M_c - M_L) Un, where M_L is the
    // lumped mass matrix, c_tau is the mass diffusion coefficient (c_tau = 1.0
    // guarantees a monotonic solution). See also the details in the function
    // header for the notation. Based on the above, the AEC, in general, is
    // computed as AEC = M_L^{-1} (M_Le - M_ce) (ctau * Un + dUh), which can be
    // obtained by subtracting the low order system from the high order system.
    // Note that the solution update is U^{n+1} = Un + dUl + lim(dUh - dUl),
    // where the last term is the limited AEC. (Think of 'lim' as the limit
    // coefficient between 0 and 1.)
    for (std::size_t j=0; j<4; ++j)
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t k=0; k<4; ++k)
          m_aec(e*4+j,c) += m[j][k] * ctau*un[c][k] / vol[N[j]];
  }

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 >
      N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
    for (std::size_t j=0; j<4; ++j) {
      // Dirichlet BCs: At nodes where Dirichlet boundary conditions (BC) are
      // set, we set the AEC to zero. This is because if the (same) BCs are
      // correctly set for both the low and the high order solution, there
      // should be no difference between the low and high order increments,
      // thus AEC = dUh - dUl = 0.
      auto b = bcdir.find(N[j]);
      if (b != end(bcdir)) {
        for (ncomp_t c=0; c<ncomp; ++c) {
          if (b->second[c].first) {
            m_aec(e*4+j,c) = 0.0;
          }
        }
      }
      // Symmetry BCs
      for (const auto& [s,nodes] : symbcnodemap) {
        auto i = nodes.find(N[j]);
        if (i != end(nodes)) {
          auto l = bnorm.find(s);
          if (l != end(bnorm)) {
            auto k = l->second.find(N[j]);
            for (const auto& vel : m_vel) {
              std::array< tk::real, 3 >
                v{ m_aec(e*4+j,vel[0]),
                   m_aec(e*4+j,vel[1]),
                   m_aec(e*4+j,vel[2]) },
                n{ k->second[0], k->second[1], k->second[2] };
              auto vn = tk::dot( v, n );
              m_aec(e*4+j,vel[0]) -= vn * n[0];
              m_aec(e*4+j,vel[1]) -= vn * n[1];
              m_aec(e*4+j,vel[2]) -= vn * n[2];
            }
          }
        }
      }
    }
  }

  // sum all positive (negative) antidiffusive element contributions to nodes
  // (Lohner: P^{+,-}_i)
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (std::size_t j=0; j<4; ++j) {
      for (ncomp_t c=0; c<ncomp; ++c) {
        P(N[j],c*2+0) += std::max( 0.0, m_aec(e*4+j,c) );
        P(N[j],c*2+1) += std::min( 0.0, m_aec(e*4+j,c) );
      }
    }
  }
}

bool
FluxCorrector::verify( std::size_t nchare,
                       const std::vector< std::size_t >& inpoel,
                       const tk::Fields& dUh,
                       const tk::Fields& dUl ) const
// *****************************************************************************
//  Verify the assembled antidiffusive element contributions (AEC)
//! \param[in] nchare Total number of host chares
//! \param[in] inpoel Mesh element connectivity
//! \param[in] dUh Increment of the high order solution
//! \param[in] dUl Increment of the low order solution
//! \return True if verification has been done
//! \details This verification only makes sense if no communication is to be
//!   done, i.e., if there is a single host chare, because the AEC assembled to
//!   mesh nodes only contains partial contributions on chare boundaries at this
//!   point. Verification in parallel would incure communication of the
//!   unlimited AEC, which in general is not necessary, so we will not do that
//!   for the sake of verification.
//! \note Client code should ensure that this function is optimized away in
//!   RELEASE mode.
// *****************************************************************************
{
  Assert( dUl.nunk() == dUh.nunk() && dUl.nprop() == dUh.nprop(),
          "Unknown array size mismatch" );

  if (nchare == 1) {
    auto ncomp = g_newinputdeck.get< newtag::ncomp >();
    tk::Fields U( dUh.nunk(), dUh.nprop() );
    U.fill( 0.0 );

    for (std::size_t e=0; e<inpoel.size()/4; ++e) {
      const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                             inpoel[e*4+2], inpoel[e*4+3] }};
      // access pointer to solution at element nodes
      std::vector< const tk::real* > u( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) u[c] = U.cptr( c );
      // scatter-add antidiffusive element contributions to nodes
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t j=0; j<4; ++j)
          U.var(u[c],N[j]) += m_aec(e*4+j,c);
    }

    // Compute maximum difference between the assembled AEC and dUh-dUl
    auto d = tk::maxdiff( U, dUh-dUl );

    // Tolerance: 10 x the linear solver tolerance for the high order solution.
    if (d.second > 1.0e-7) {
      const auto& duh = dUh.vec();
      const auto& dul = dUl.vec();
      const auto& u = U.vec();
      std::stringstream ss;
      ss << "maximum difference at mesh node " << d.first << ": " << d.second
         << ", dUh:" << duh[d.first] << ", dUl:" << dul[d.first]
         << ", dUh-dUl:" << duh[d.first] - dul[d.first]
         << ", AEC:" << u[d.first];
      Throw( "Assembled AEC does not equal dUh-dUl: " + ss.str() );
    }

    return true;
  }

  return false;
}

tk::Fields
FluxCorrector::diff( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const tk::Fields& Un ) const
// *****************************************************************************
//  Compute mass diffusion contribution to the RHS of the low order system
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] Un Solution at the previous time step
//! \return Mass diffusion contribution to the RHS of the low order system
// *****************************************************************************
{
  auto ncomp = g_newinputdeck.get< newtag::ncomp >();
  auto ctau = g_newinputdeck.get< newtag::ctau >();

  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  tk::Fields D( Un.nunk(), Un.nprop() );
  D.fill( 0.0 );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
       N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
     // compute element Jacobi determinant
     const std::array< tk::real, 3 >
       ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
       ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
       da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
     const auto J = tk::triple( ba, ca, da );   // J = 6V
     Assert( J > 0, "Element Jacobian non-positive" );

     // lumped - consistent mass
     std::array< std::array< tk::real, 4 >, 4 > m;       // nnode*nnode [4][4]
     m[0][0] = m[1][1] = m[2][2] = m[3][3] = 3.0*J/120.0;// diagonal
     m[0][1] = m[0][2] = m[0][3] =                       // off-diagonal
     m[1][0] = m[1][2] = m[1][3] =
     m[2][0] = m[2][1] = m[2][3] =
     m[3][0] = m[3][1] = m[3][2] = -J/120.0;

     // access solution at element nodes at time n
     std::vector< std::array< tk::real, 4 > > un( ncomp );
     for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, N );
     // access pointer to mass diffusion right hand side at element nodes
     std::vector< const tk::real* > d( ncomp );
     for (ncomp_t c=0; c<ncomp; ++c) d[c] = D.cptr( c );

     // scatter-add mass diffusion element contributions to rhs nodes
     for (std::size_t a=0; a<4; ++a) {
       for (ncomp_t c=0; c<ncomp; ++c)
         for (std::size_t b=0; b<4; ++b)
           D.var(d[c],N[a]) -= ctau * m[a][b] * un[c][b];
     }
  }

  return D;
}

void
FluxCorrector::alw( const std::vector< std::size_t >& inpoel,
                    const tk::Fields& Un,
                    const tk::Fields& Ul,
                    tk::Fields& Q ) const
// *****************************************************************************
//  Compute the maximum and minimum unknowns of elements surrounding nodes
//! \param[in] inpoel Mesh element connectivity
//! \param[in] Un Solution at the previous time step
//! \param[in] Ul Low order solution
//! \param[in,out] Q Maximum and mimimum unknowns of elements surrounding nodes
// *****************************************************************************
{
  Assert( Q.nunk() == Un.nunk() && Q.nprop() == Un.nprop()*2, "Max and min "
          "unknowns of elements surrounding nodes array size mismatch" );

  auto ncomp = g_newinputdeck.get< newtag::ncomp >();
  auto clip = g_newinputdeck.get< newtag::fctclip >();

  // compute maximum and minimum nodal values of all elements (Lohner: u^*_el)
  tk::Fields S( inpoel.size()/4, ncomp*2 );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      S(e,c*2+0) = -std::numeric_limits< tk::real >::max();
      S(e,c*2+1) = std::numeric_limits< tk::real >::max();
      for (std::size_t j=0; j<4; ++j) {
        // compute maximum and minimum nodal values of Ul and Un (Lohner: u^*_i)
        auto jmax = clip ? Ul(N[j],c) : std::max(Ul(N[j],c), Un(N[j],c));
        auto jmin = clip ? Ul(N[j],c) : std::min(Ul(N[j],c), Un(N[j],c));
        if (jmax > S(e,c*2+0)) S(e,c*2+0) = jmax;
        if (jmin < S(e,c*2+1)) S(e,c*2+1) = jmin;
      }
    }
  }

  // compute maximum and mimimum unknowns of all elements surrounding each node
  // (Lohner: u^{max,min}_i)
  const auto esup = tk::genEsup( inpoel, 4 );
  for (std::size_t p=0; p<Un.nunk(); ++p) {
    for (auto e : tk::Around(esup,p)) {
      for (ncomp_t c=0; c<ncomp; ++c) {
        if (S(e,c*2+0) > Q(p,c*2+0)) Q(p,c*2+0) = S(e,c*2+0);
        if (S(e,c*2+1) < Q(p,c*2+1)) Q(p,c*2+1) = S(e,c*2+1);
      }
    }
  }
}

void
FluxCorrector::lim( const std::vector< std::size_t >& inpoel,
                    const std::unordered_map< std::size_t,
                            std::vector< std::pair< bool, tk::real > > >& bcdir,
                    const tk::Fields& P,
                    const tk::Fields& Ul,
                    tk::Fields& Q,
                    tk::Fields& A ) const
// *****************************************************************************
// Compute limited antiffusive element contributions and apply to mesh nodes
//! \param[in] inpoel Mesh element connectivity
//! \param[in] bcdir Vector of pairs of bool and boundary condition value
//!   associated to mesh node IDs at which to set Dirichlet boundary conditions.
//! \param[in] P The sums of all positive (negative) AECs to nodes
//! \param[in] Ul Low order solution
//! \param[in,out] Q The maximum and mimimum unknowns of elements surrounding
//!   each node
//! \param[in,out] A Limited antidiffusive element contributions scatter-added
//!   to nodes
//! \note Q is also overwritten to avoid using temporary memory
// *****************************************************************************
{
  Assert( P.nunk() == Q.nunk() && P.nprop() == Q.nprop(), "Size mismatch" );
  Assert( P.nunk() == Ul.nunk() && P.nprop() == Ul.nprop()*2, "Size mismatch" );
  Assert( A.nunk() == Ul.nunk() && A.nprop() == Ul.nprop(), "Size mismatch" );

  auto ncomp = g_newinputdeck.get< newtag::ncomp >();

  // compute the maximum and minimum increments and decrements nodal solution
  // values are allowed to achieve (Lohner: Q^{+,-}_i)
  for (std::size_t p=0; p<Ul.nunk(); ++p)
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c*2+0) -= Ul(p,c);
      Q(p,c*2+1) -= Ul(p,c);
    }

  auto eps = g_newinputdeck.get< newtag::fcteps >();

  // compute the ratios of positive and negative element contributions that
  // ensure monotonicity (Lohner: R^{+,-})
  for (std::size_t p=0; p<P.nunk(); ++p) {
    for (ncomp_t c=0; c<ncomp; ++c) {

      if (P(p,c*2+0) < eps)
        Q(p,c*2+0) = 1.0;
      else
        Q(p,c*2+0) = std::min(1.0,Q(p,c*2+0)/P(p,c*2+0));

      if (P(p,c*2+1) > -eps)
        Q(p,c*2+1) = 1.0;
      else
        Q(p,c*2+1) = std::min(1.0,Q(p,c*2+1)/P(p,c*2+1));

    }
  }

  // calculate limit coefficient for all elements (Lohner: C_el)
  tk::Fields C( inpoel.size()/4, ncomp );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      std::array< tk::real, 4 > R;
      for (std::size_t j=0; j<4; ++j) {

        if (std::abs(m_aec(e*4+j,c)) < eps)
          R[j] = 1.0;
        else if (m_aec(e*4+j,c) > 0.0)
          R[j] = Q(N[j],c*2+0);
        else
          R[j] = Q(N[j],c*2+1);

      }
      C(e,c) = *std::min_element( begin(R), end(R) );
      // if all vertices happened to be on a Dirichlet boundary, ignore limiting
      if (C(e,c) > 1.0) C(e,c) = 1.0;
      Assert( C(e,c) > -eps && C(e,c) < 1.0+eps,
              "0 <= AEC <= 1.0 failed: C = " + std::to_string(C(e,c)) );
    }
  }

  // System limiting
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    tk::real cs = 1.0;
    for (auto i : m_sys) if (C(e,i) < cs) cs = C(e,i);
    for (auto i : m_sys) C(e,i) = cs;
  }

  // save the limited antidiffusive element contributions (Lohner: AEC^c)
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};

    // access pointer to solution at element nodes
    std::vector< const tk::real* > a( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) a[c] = A.cptr( c );

    // Scatter-add limited antidiffusive element contributions to nodes. At
    // nodes where Dirichlet boundary conditions are set, the AECs are set to
    // zero so the limit coefficient has no effect. This yields no increment for
    // those nodes. See the detailed discussion when computing the AECs.
    for (std::size_t j=0; j<4; ++j) {
      auto b = bcdir.find( N[j] );    // Dirichlet BC
      for (ncomp_t c=0; c<ncomp; ++c) {
        if (b != end(bcdir) && b->second[c].first) {
          A.var(a[c],N[j]) += m_aec(e*4+j,c);
        } else {
          A.var(a[c],N[j]) += C(e,c) * m_aec(e*4+j,c);
        }
      }
    }
  }
}

std::tuple< std::vector< std::string >,
            std::vector< std::vector< tk::real > > >
FluxCorrector::fields( const std::vector< std::size_t >& /*inpoel*/ ) const
// *****************************************************************************
//  Collect mesh output fields from FCT
//! \return Names and fields in mesh cells
// *****************************************************************************
{
  using tuple_t = std::tuple< std::vector< std::string >,
                              std::vector< std::vector< tk::real > > >;
  return tuple_t{};
}
