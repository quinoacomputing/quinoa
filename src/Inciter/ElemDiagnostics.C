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
                          const tk::Fields& u )
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization proxy to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] u Current solution vector
//! \return True if diagnostics have been computed
// *****************************************************************************
{
  // Optionally collect diagnostics and send for aggregation across all workers

  // Query after how many time steps user wants to dump diagnostics
  auto diagfreq = g_inputdeck.get< tag::interval, tag::diag >();

  if ( !((d.It()+1) % diagfreq) ) {     // if remainder, don't dump

    // Collect analytical solutions (if available) from all PDEs. Note that
    // calling the polymorphic PDE::initialize() is assumed to evaluate the
    // analytical solution for a PDE. For those PDE problems that have
    // analytical solutions, this is the same as used for setting the initial
    // conditions, since if the analytical solution is available for a problem,
    // it is (so far anyway) always initialized to its analytical solution and
    // that is done by calling PDE::initialize(). If the analytical solution is
    // a function of time, that is already incorporated in setting initial
    // conditions. For those PDEs where the analytical solution is not
    // available, initialize() returns the initial conditions (obviously), and
    // thus the "error", defined between this "analytical" and numerical
    // solution will be a measure of the "distance" between the initial
    // condition and the current numerical solution. This is not necessarily
    // useful, but simplies the logic because all PDEs can be treated as being
    // able to compute an error based on some "analytical" solution, which is
    // really the initial condition.
    auto a = u;
    for (const auto& eq : g_dgpde)
      eq.initialize( geoElem, a, d.T()+d.Dt() );

    // Prepare for computing diagnostics. Diagnostics are defined as some norm,
    // .e.g., L2 norm, of a quantity, computed in mesh nodes, A, as ||A||_2 =
    // sqrt[ sum_i(A_i)^2 V_i ], where the sum is taken over all mesh nodes and
    // V_i is the nodal volume. We send multiple sets of quantities to the host
    // for aggregation across the whole mesh: (1) the numerical solutions of all
    // components of all PDEs, and their error, defined as A_i = (a_i - n_i),
    // where a_i and n_i are the analytical and numerical solutions at node i,
    // respectively. The final aggregated solution will end up in
    // Transporter::diagnostics().

    // Diagnostics vector (of vectors) during aggregation:
    // 0: L2-norm of all scalar components of the numerical solution
    // 1: L2-norm of all scalar components of the numerical-analytic solution
    // 2: Linf-norm of all scalar components of the numerical-analytic solution
    auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop()/ndof, 0.0 ) );

    // Put in norms sweeping our mesh chunk
    for (std::size_t i=0; i<u.nunk()-nchGhost; ++i) {
      // Compute sum for L2 norm of the numerical solution
      for (std::size_t c=0; c<u.nprop()/ndof; ++c)
      {
        auto mark = c*ndof;
        diag[L2SOL][c] += u(i,mark,0) * u(i,mark,0) * geoElem(i,0,0);
      }
      // Compute sum for L2 norm of the numerical-analytic solution
      for (std::size_t c=0; c<u.nprop()/ndof; ++c)
      {
        auto mark = c*ndof;
        diag[L2ERR][c] +=
          (u(i,mark,0)-a(i,mark,0)) * (u(i,mark,0)-a(i,mark,0)) * geoElem(i,0,0);
      }
      // Compute max for Linf norm of the numerical-analytic solution
      for (std::size_t c=0; c<u.nprop()/ndof; ++c) {
        auto mark = c*ndof;
        auto err = std::abs( u(i,mark,0) - a(i,mark,0) );
        if (err > diag[LINFERR][c]) diag[LINFERR][c] = err;
      }
    }

    // Append diagnostics vector with metadata on the current time step
    // 3: Current iteration count (only the first entry is used)
    // 4: Current physical time (only the first entry is used)
    // 5: Current physical time step size (only the first entry is used)
    diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
    diag[TIME][0] = d.T() + d.Dt();
    diag[DT][0] = d.Dt();

    // Contribute to diagnostics
    auto stream = serialize( diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;
  }

  return false;
}

bool
ElemDiagnostics::computep1( Discretization& d,
                            const std::size_t nchGhost,
                            const tk::Fields& geoElem,
                            const std::vector< std::size_t >& inpoel,
                            const tk::UnsMesh::Coords& coord,
                            const tk::Fields& u )
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc. for dgp1
//! \param[in] d Discretization proxy to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] u Current solution vector
//! \return True if diagnostics have been computed
// *****************************************************************************
{
  // Optionally collect diagnostics and send for aggregation across all workers

  // Query after how many time steps user wants to dump diagnostics
  auto diagfreq = g_inputdeck.get< tag::interval, tag::diag >();

  if ( !((d.It()+1) % diagfreq) ) {     // if remainder, don't dump

    // Prepare for computing diagnostics. Diagnostics are defined as some norm,
    // .e.g., L2 norm, of a quantity, computed in mesh nodes, A, as ||A||_2 =
    // sqrt[ sum_i(A_i)^2 V_i ], where the sum is taken over all mesh nodes and
    // V_i is the nodal volume. We send multiple sets of quantities to the host
    // for aggregation across the whole mesh: (1) the numerical solutions of all
    // components of all PDEs, and their error, defined as A_i = (a_i - n_i),
    // where a_i and n_i are the analytical and numerical solutions at node i,
    // respectively. The final aggregated solution will end up in
    // Transporter::diagnostics().

    // Diagnostics vector (of vectors) during aggregation:
    // 0: L2-norm of all scalar components of the numerical solution
    // 1: L2-norm of all scalar components of the numerical-analytic solution
    // 2: Linf-norm of all scalar components of the numerical-analytic solution
    auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop()/ndof, 0.0 ) );

    // Quadrature points
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize(5,0);
    coordgp[1].resize(5,0);
    coordgp[2].resize(5,0);

    wgp.resize(5,0);

    coordgp[0][0] = 0.25;
    coordgp[1][0] = 0.25;
    coordgp[2][0] = 0.25;
    wgp[0]        = -12.0/15.0;

    coordgp[0][1] = 1.0/6.0;
    coordgp[1][1] = 1.0/6.0;
    coordgp[2][1] = 1.0/6.0;
    wgp[1]        = 9.0/20.0;

    coordgp[0][2] = 0.5;
    coordgp[1][2] = 1.0/6.0;
    coordgp[2][2] = 1.0/6.0;
    wgp[2]        = 9.0/20.0;

    coordgp[0][3] = 1.0/6.0;
    coordgp[1][3] = 0.5;
    coordgp[2][3] = 1.0/6.0;
    wgp[3]        = 9.0/20.0;

    coordgp[0][4] = 1.0/6.0;
    coordgp[1][4] = 1.0/6.0;
    coordgp[2][4] = 0.5;
    wgp[4]        = 9.0/20.0;

    // Put in norms sweeping our mesh chunk
    for (std::size_t i=0; i<u.nunk()-nchGhost; ++i) {

      auto vole = geoElem(i,0,0);

      auto x1 = coord[0][ inpoel[4*i]   ];
      auto y1 = coord[1][ inpoel[4*i]   ];
      auto z1 = coord[2][ inpoel[4*i]   ];

      auto x2 = coord[0][ inpoel[4*i+1] ];
      auto y2 = coord[1][ inpoel[4*i+1] ];
      auto z2 = coord[2][ inpoel[4*i+1] ];

      auto x3 = coord[0][ inpoel[4*i+2] ];
      auto y3 = coord[1][ inpoel[4*i+2] ];
      auto z3 = coord[2][ inpoel[4*i+2] ];

      auto x4 = coord[0][ inpoel[4*i+3] ];
      auto y4 = coord[1][ inpoel[4*i+3] ];
      auto z4 = coord[2][ inpoel[4*i+3] ];

      // Gaussian quadrature
      for (std::size_t igp=0; igp<5; ++igp)
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
          s = eq.analyticalSol( xgp, ygp, zgp, d.T()+d.Dt() );

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

    // Append diagnostics vector with metadata on the current time step
    // 3: Current iteration count (only the first entry is used)
    // 4: Current physical time (only the first entry is used)
    // 5: Current physical time step size (only the first entry is used)
    diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
    diag[TIME][0] = d.T() + d.Dt();
    diag[DT][0] = d.Dt();

    // Contribute to diagnostics
    auto stream = serialize( diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;
  }

  return false;
}
