// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "Macro.h"
#include "Vector.h"
#include "DerivedData.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

using inciter::FluxCorrector;

void
FluxCorrector::aec( const std::array< std::vector< tk::real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const std::vector< tk::real >& vol,
                    const std::unordered_map< std::size_t,
                            std::vector< std::pair< bool, tk::real > > >& bc,
                    const std::vector< std::size_t >& gid,
                    const tk::Fields& dUh,
                    const tk::Fields& Un,
                    tk::Fields& P )
// *****************************************************************************
//  Compute antidiffusive element contributions (AEC)
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] vol Volume associated to mesh nodes
//! \param[in] bc Vector of pairs of bool and boundary condition value
//!   associated to mesh node IDs at which to set Dirichlet boundary conditions.
//!   Note that this BC data structure must include boundary conditions set
//!   across all PEs, not just the ones need to be set on this PE.
//! \param[in] gid Local to global node ID mapping
//! \param[in] dUh Increment of the high order solution
//! \param[in] Un Solution at the previous time step stage
//! \param[in,out] P The sums of positive (negative) AECs to nodes
//! \details The antidiffusive element contributions (AEC) are defined as the
//!   difference between the high and low order solution, where the high order
//!   solution is obtained from consistent mass Taylor-Galerkin discretization
//!   and the low order solution is lumped mass Taylor-Galerkin + diffusion.
//!   Note that AEC is not directly computed as dUh - dUl (although that could
//!   also be done), but as AEC =  M_L^{-1} (M_Le - M_ce) (ctau * Un + dUh),
//!   where
//!    * M_ce is the element's consistent mass matrix,
//!    * M_Le is the element's lumped mass matrix,
//!    * ctau is the mass diffusion coefficient on the rhs of the low order
//!      solution, see also FluxCorrector::diff(),
//!    * Un is the solution at the previous time step stage,
//!    * dUh is the increment of the high order solution, and
//!    * M_L^{-1} is the inverse of the assembled lumped mass matrix, i.e., the
//!      volume associated to a mesh node by summing the quarter of the element
//!      volumes surrounding the node. Note that this is the correct node volume
//!      taking into account that some nodes are on chare boundaries. See also
//!      Carrier::vol().
//! \see Löhner, R., Morgan, K., Peraire, J. and Vahdati, M. (1987), Finite
//!   element flux-corrected transport (FEM–FCT) for the Euler and Navier–Stokes
//!   equations. Int. J. Numer. Meth. Fluids, 7: 1093–1109.
//!   doi:10.1002/fld.1650071007
//! \author J. Bakosi
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

  Assert( vol.size() == coord[0].size(), "Nodal volume vector size mismatch" );
  Assert( m_aec.nunk() == inpoel.size() && m_aec.nprop() == ncomp,
          "AEC and mesh connectivity size mismatch" );
  Assert( Un.nunk() == dUh.nunk() && Un.nprop() == dUh.nprop(),
          "Unknown array size mismatch" );
  Assert( P.nunk() == dUh.nunk() && P.nprop() == dUh.nprop()*2,
          "Sums of positive (negative) AECs to nodes array size mismatch" );

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
    for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, 0, N );
    // access high-order solution increment at element nodes
    std::vector< std::array< tk::real, 4 > > duh( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) duh[c] = dUh.extract( c, 0, N );

    // Compute antidiffusive element contributions (AEC). The high order system
    // is M_c * dUh = r, where M_c is the consistent mass matrix and r is the
    // high order right hand side. The low order system is constructed from the
    // high order one by lumping the consistent mass matrix and adding mass
    // diffusion: M_L * dUl = r + c_tau * (M_c - M_L) Un, where M_L is the
    // lumped mass matrix, c_tau is the mass diffusion coefficient (c_tau = 1.0
    // guarantees a monotonic solution). See also the details in the function
    // header for the notaion. Based on the above, the AEC, in general, is
    // computed as AEC = M_L^{-1} (M_Le - M_ce) (ctau * Un + dUh), which can be
    // obtained by subtracting the low order system from the high order system.
    // Note that the solution update is U^{n+1} = Un + dUl + lim(dUh - dUl),
    // where the last terms is the limited AEC. (Think of 'lim' as the limit
    // coefficient between 0 and 1.)
    for (std::size_t j=0; j<4; ++j)
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t k=0; k<4; ++k)
          m_aec(e*4+j,c,0) += m[j][k] * (ctau*un[c][k] + duh[c][k]) / vol[N[j]];
  }

  // At nodes where Dirichlet boundary conditions (BC) are set, we set the AEC
  // to zero. Since the right hand side of the low order solution is also set to
  // zero at nodes where Dirichlet BCs are set, this properly enforces no
  // increment at BC nodes. See also LinSysMerger::auxsolve().
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (std::size_t j=0; j<4; ++j) {
      auto b = bc.find( gid[ N[j] ] );
      if (b != end(bc))
        for (ncomp_t c=0; c<ncomp; ++c)
          if (b->second[c].first)
            m_aec(e*4+j,c,0) = 0.0;
    }
  }

  // sum all positive (negative) antidiffusive element contributions to nodes
  // (Lohner: P^{+,-}_i)
  const auto esup = tk::genEsup( inpoel, 4 );
  for (std::size_t p=0; p<dUh.nunk(); ++p)
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i) {
       const auto e = esup.first[i];
       // decide which node's contribution we need (A, B, C, or D)
       std::size_t n = 3;
            if (inpoel[e*4+0] == p) n = 0;
       else if (inpoel[e*4+1] == p) n = 1;
       else if (inpoel[e*4+2] == p) n = 2;
       // add up element contributions to node p
       for (ncomp_t c=0; c<ncomp; ++c) {
         P(p,c*2+0,0) += std::max( 0.0, m_aec(e*4+n,c,0) );
         P(p,c*2+1,0) += std::min( 0.0, m_aec(e*4+n,c,0) );
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
//! \note This function is optimized away in RELEASE mode, see carrier.ci and
//!   Carrier::verify().
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( dUl.nunk() == dUh.nunk() && dUl.nprop() == dUh.nprop(),
          "Unknown array size mismatch" );

  if (nchare == 1) {
    auto ncomp = g_inputdeck.get< tag::component >().nprop();
    tk::Fields U( dUh.nunk(), dUh.nprop() );
    U.fill( 0.0 );

    for (std::size_t e=0; e<inpoel.size()/4; ++e) {
      const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                             inpoel[e*4+2], inpoel[e*4+3] }};
      // access pointer to solution at element nodes
      std::vector< const tk::real* > u( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) u[c] = U.cptr( c, 0 );
      // scatter-add antidiffusive element contributions to nodes
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t j=0; j<4; ++j)
          U.var(u[c],N[j]) += m_aec(e*4+j,c,0);
    }

    // Compute maximum difference between the assembled AEC and dUh-dUl
    auto d = tk::maxdiff( U, dUh-dUl );

    // Tolerance: 10 x the linear solver tolerance for the high order solution.
    if (d.second > 1.0e-7) {
      const auto& duh = dUh.data();
      const auto& dul = dUl.data();
      const auto& u = U.data();
      std::stringstream ss;
      ss << "maximum difference at " << d.first << ": " << d.second
         << ", dUh:" << duh[d.first] << ", dUl:" << dul[d.first]
         << ", dUh-dUl:" << duh[d.first] - dul[d.first]
         << ", AEC:" << u[d.first] << '\n';
      Throw( "Assembled AEC does not equal dUh-dUl: " + ss.str() );
    }

    return true;
  }

  return false;
}

tk::Fields
FluxCorrector::lump( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel ) const
// *****************************************************************************
//  Compute lumped mass matrix left hand side for the low order system
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \return Lumped mass matrix
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  tk::Fields L( coord[0].size(), ncomp );
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
    Assert( J > 0, "Element Jacobian non-positive" );

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

tk::Fields
FluxCorrector::diff( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const tk::Fields& Un ) const
// *****************************************************************************
//  Compute mass diffusion contribution to the RHS of the low order system
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] Un Solution at the previous time step stage
//! \return Mass diffusion contribution to the RHS of the low order system
//! \author J. Bakosi
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  auto ctau = g_inputdeck.get< tag::discr, tag::ctau >();

  tk::Fields D( Un.nunk(), Un.nprop() );
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
    for (ncomp_t c=0; c<ncomp; ++c) un[c] = Un.extract( c, 0, N );
    // access pointer to mass diffusion right hand side at element nodes
    std::vector< const tk::real* > d( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) d[c] = D.cptr( c, 0 );

    // scatter-add mass diffusion element contributions to rhs nodes
    for (std::size_t j=0; j<4; ++j)
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t k=0; k<4; ++k)
          D.var(d[c],N[j]) -= ctau * m[j][k] * un[c][k];
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
//! \param[in] Un Solution at the previous time step stage
//! \param[in] Ul Low order solution
//! \param[in,out] Q Maximum and mimimum unknowns of elements surrounding nodes
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( Q.nunk() == Un.nunk() && Q.nprop() == Un.nprop()*2, "Max and min "
          "unknowns of elements surrounding nodes array size mismatch" );

  // compute maximum and minimum nodal values of Ul and Un (Lohner: u^*_i)
  auto Smax = tk::max( Ul, Un );
  auto Smin = tk::min( Ul, Un );

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  // compute maximum and minimum nodal values of all elements (Lohner: u^*_el)
  tk::Fields S( inpoel.size()/4, ncomp*2 );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      S(e,c*2+0,0) = std::max( Smax(N[3],c,0),
                       std::max( Smax(N[2],c,0),
                         std::max( Smax(N[0],c,0), Smax(N[1],c,0) ) ) );
      S(e,c*2+1,0) = std::min( Smin(N[3],c,0),
                       std::min( Smin(N[2],c,0),
                         std::min( Smin(N[0],c,0), Smin(N[1],c,0) ) ) );
    }
  }

  // compute maximum and mimimum unknowns of all elements surrounding each node
  // (Lohner: u^{max,min}_i)
  const auto esup = tk::genEsup( inpoel, 4 );
  for (std::size_t p=0; p<Un.nunk(); ++p)
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i) {
      const auto e = esup.first[i];
      for (ncomp_t c=0; c<ncomp; ++c) {
        if (S(e,c*2+0,0) > Q(p,c*2+0,0)) Q(p,c*2+0,0) = S(e,c*2+0,0);
        if (S(e,c*2+1,0) < Q(p,c*2+1,0)) Q(p,c*2+1,0) = S(e,c*2+1,0);
      }
    }
}

void
FluxCorrector::lim( const std::vector< std::size_t >& inpoel,
                    const tk::Fields& P,
                    const tk::Fields& Ul,
                    tk::Fields& Q,
                    tk::Fields& A ) const
// *****************************************************************************
// Compute limited antiffusive element contributions and apply to mesh nodes
//! \param[in] inpoel Mesh element connectivity
//! \param[in] P The sums of all positive (negative) AECs to nodes
//! \param[in] Ul Low order solution
//! \param[inout] Q The maximum and mimimum unknowns of elements surrounding
//!   each node
//! \param[in,out] Limited antidiffusive element contributions scatter-added to
//!   nodes
//! \note Q is also overwritten to avoid using temporary memory
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( P.nunk() == Q.nunk() && P.nprop() == Q.nprop(), "Size mismatch" );
  Assert( P.nunk() == Ul.nunk() && P.nprop() == Ul.nprop()*2, "Size mismatch" );
  Assert( A.nunk() == Ul.nunk() && A.nprop() == Ul.nprop(), "Size mismatch" );

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  // compute the maximum and minimum increments and decrements nodal solution
  // values are allowed to achieve (Lohner: Q^{+,-}_i)
  for (std::size_t p=0; p<Ul.nunk(); ++p)
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c*2+0,0) -= Ul(p,c,0);
      Q(p,c*2+1,0) -= Ul(p,c,0);
    }

  // compute the ratios of positive and negative element contributions that
  // ensure monotonicity (Lohner: R^{+,-})
  for (std::size_t p=0; p<P.nunk(); ++p)
    for (ncomp_t c=0; c<ncomp; ++c) {
      Q(p,c*2+0,0) =
        (P(p,c*2+0,0) > 0.0 ? std::min(1.0,Q(p,c*2+0,0)/P(p,c*2+0,0)) : 0.0);
      Q(p,c*2+1,0) =
        (P(p,c*2+1,0) < 0.0 ? std::min(1.0,Q(p,c*2+1,0)/P(p,c*2+1,0)) : 0.0);
    }

  // calculate limit coefficient for all elements (Lohner: C_el)
  tk::Fields C( inpoel.size()/4, ncomp );
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    for (ncomp_t c=0; c<ncomp; ++c) {
      std::array< tk::real, 4 > R;
      for (std::size_t j=0; j<4; ++j)
        R[j] = m_aec(e*4+j,c,0) > 0.0 ? Q(N[j],c*2+0,0) : Q(N[j],c*2+1,0);
      C(e,c,0) = *std::min_element( begin(R), end(R) );
      Assert( C(e,c,0) > -std::numeric_limits< tk::real >::epsilon() &&
              C(e,c,0) < 1.0+std::numeric_limits< tk::real >::epsilon(),
              "0 <= AEC <= 1.0 failed" );
    }
  }

  // save the limited antidiffusive element contributions (Lohner: AEC^c)
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};

    // access pointer to solution at element nodes
    std::vector< const tk::real* > a( ncomp );
    for (ncomp_t c=0; c<ncomp; ++c) a[c] = A.cptr( c, 0 );

    // Scatter-add limited antidiffusive element contributions to nodes. At
    // nodes where Dirichlet boundary conditions are set, we ignore the limit
    // coefficient. This yields no increment for those nodes. See the detailed
    // discussion when computing the AECs.
    for (std::size_t j=0; j<4; ++j)
      for (ncomp_t c=0; c<ncomp; ++c)
        A.var(a[c],N[j]) += C(e,c,0) * m_aec(e*4+j,c,0);
  }
}
