// *****************************************************************************
/*!
  \file      src/Inciter/ElemDiagnostics.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     ElemDiagnostics class for collecting element diagnostics
  \details   ElemDiagnostics class for collecting element diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************

#include "DGPDE.h"
#include "ElemDiagnostics.h"
#include "DiagReducer.h"
#include "Discretization.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DGPDE > g_dgpde;

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::ElemDiagnostics;

void
ElemDiagnostics::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details This routine is supposed to be called from a Charm++ nodeinit
//!   routine. Since the runtime system executes nodeinit routines exactly once
//!   on every logical node early on in the Charm++ init sequence, they must be
//!   static as they are called without an object. See also: Section
//!   "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  DiagMerger = CkReduction::addReducer( mergeDiag );
}

bool
ElemDiagnostics::compute( Discretization& d,
                          const std::size_t nchGhost,
                          const tk::Fields& geoElem,
                          const tk::Fields& u ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization base class to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] u Current solution vector
//! \return True if diagnostics have been computed
//! \details Diagnostics are defined as some norm, e.g., L2 norm, of a quantity,
//!    computed in mesh elements, A, as ||A||_2 = sqrt[ sum_i(A_i)^2 V_i ],
//!    where the sum is taken over all mesh elements and V_i is the cell volume.
//!    We send multiple sets of quantities to the host for aggregation across
//!    the whole mesh. The final aggregated solution will end up in
//!    Transporter::diagnostics(). Aggregation of the partially computed
//!    diagnostics is done via potentially different policies for each field.
//! \see inciter::mergeDiag(), src/Inciter/Diagnostics.h
// *****************************************************************************
{
  // Optionally collect diagnostics and send for aggregation across all workers

  // Query after how many time steps user wants to dump diagnostics
  auto diagfreq = g_inputdeck.get< tag::interval, tag::diag >();

  if ( !((d.It()+1) % diagfreq) ) {  // if remainder, don't compute diagnostics

    // Query number of degrees of freedom from user's setting
    const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();

    // Diagnostics vector (of vectors) during aggregation. See
    // Inciter/Diagnostics.h.
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop()/ndof, 0.0 ) );
  
    if (ndof == 1)
      computeP0( d, nchGhost, geoElem, u, diag );
    else if (ndof == 4)
      computeP1( d, nchGhost, geoElem, u, diag );
    else Throw( "ElemDiagnostics::compute() not defined for NDOF=" +
                 std::to_string(ndof) );

    // Append diagnostics vector with metadata on the current time step
    // ITER: Current iteration count (only the first entry is used)
    // TIME: Current physical time (only the first entry is used)
    // DT: Current physical time step size (only the first entry is used)
    diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
    diag[TIME][0] = d.T() + d.Dt();
    diag[DT][0] = d.Dt();

    // Contribute to diagnostics
    auto stream = serialize( diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;        // diagnostics have been computed

  }

  return false;         // diagnostics have not been computed
}

void
ElemDiagnostics::computeP0( Discretization& d,
                            const std::size_t nchGhost,
                            const tk::Fields& geoElem,
                            const tk::Fields& u,
                            std::vector< std::vector< tk::real > >& diag ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc. for DG(P0)
//! \param[in] d Discretization base class to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] u Current solution vector
//! \param[in,out] diag Diagnostics vector
// *****************************************************************************
{
  for (std::size_t i=0; i<u.nunk()-nchGhost; ++i) {

    // Compute sum for L2 norm of the numerical solution
    for (std::size_t c=0; c<u.nprop(); ++c)
      diag[L2SOL][c] += u(i,c,0) * u(i,c,0) * geoElem(i,0,0);

    // Query and collect analytic solution for all components of all PDEs
    // integrated at cell centroids
    std::vector< tk::real > a;
    const auto xc = geoElem(i,1,0);
    const auto yc = geoElem(i,2,0);
    const auto zc = geoElem(i,3,0);
    for (const auto& eq : g_dgpde) {
      auto s = eq.analyticSolution( xc, yc, zc, d.T()+d.Dt() );
      std::move( begin(s), end(s), std::back_inserter(a) );
    }
    Assert( a.size() == u.nprop(), "Size mismatch" );

    // Compute sum for L2 norm of the numerical-analytic solution
    for (std::size_t c=0; c<u.nprop(); ++c)
      diag[L2ERR][c] += (u(i,c,0)-a[c]) * (u(i,c,0)-a[c]) * geoElem(i,0,0);
    // Compute max for Linf norm of the numerical-analytic solution
    for (std::size_t c=0; c<u.nprop(); ++c) {
      auto err = std::abs( u(i,c,0) - a[c] );
      if (err > diag[LINFERR][c]) diag[LINFERR][c] = err;
    }

  }
}

void
ElemDiagnostics::computeP1( Discretization& d,
                            const std::size_t nchGhost,
                            const tk::Fields& geoElem,
                            const tk::Fields& u,
                            std::vector< std::vector< tk::real > >& diag ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc. for DG(P1)
//! \param[in] d Discretization base class to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] u Current solution vector
//! \return True if diagnostics have been computed
// *****************************************************************************
{
  const auto& inpoel = d.Inpoel();
  const auto& coord = d.Coord();

  // Number of integration points
  constexpr std::size_t NG = 4;

  // Quadrature points
  std::array< std::array< tk::real, NG >, 3 > coordgp;
  std::array< tk::real, NG > wgp;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  GaussQuadratureTet( coordgp, wgp );

  // Query number of degrees of freedom from user's setting
  const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();

  // Put in norms sweeping our mesh chunk
  for (std::size_t i=0; i<u.nunk()-nchGhost; ++i) {

    auto vole = geoElem(i,0,0);

    auto x1 = cx[ inpoel[4*i]   ];
    auto y1 = cy[ inpoel[4*i]   ];
    auto z1 = cz[ inpoel[4*i]   ];

    auto x2 = cx[ inpoel[4*i+1] ];
    auto y2 = cy[ inpoel[4*i+1] ];
    auto z2 = cz[ inpoel[4*i+1] ];

    auto x3 = cx[ inpoel[4*i+2] ];
    auto y3 = cy[ inpoel[4*i+2] ];
    auto z3 = cz[ inpoel[4*i+2] ];

    auto x4 = cx[ inpoel[4*i+3] ];
    auto y4 = cy[ inpoel[4*i+3] ];
    auto z4 = cz[ inpoel[4*i+3] ];

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
      auto B2 = 2.0 * coordgp[0][igp] + coordgp[1][igp] + coordgp[2][igp]
                - 1.0;
      auto B3 = 3.0 * coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B4 = 4.0 * coordgp[2][igp] - 1.0;

      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      auto xgp = x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4;

      auto wt = vole * wgp[igp];

      std::vector< tk::real > s;

      for (const auto& eq : g_dgpde)
        s = eq.analyticSolution( xgp, ygp, zgp, d.T()+d.Dt() );

      for (std::size_t c=0; c<u.nprop()/ndof; ++c)
      {
        auto mark = c*ndof;
        auto ugp =   u(i, mark,   0)
                   + u(i, mark+1, 0) * B2
                   + u(i, mark+2, 0) * B3
                   + u(i, mark+3, 0) * B4;

        // Compute sum for L2 norm of the numerical solution
        diag[L2SOL][c] += wt * ugp * ugp;

        // Compute sum for L2 norm of the numerical-analytic solution
        diag[L2ERR][c] += wt * (ugp-s[c]) * (ugp-s[c]);

        // Compute max for Linf norm of the numerical-analytic solution
        auto err = std::abs( ugp - s[c] );
        if (err > diag[LINFERR][c]) diag[LINFERR][c] = err;
      }
    }
  }
}
