// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     NodeDiagnostics class for collecting nodal diagnostics
  \details   NodeDiagnostics class for collecting nodal diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************

#include "CGPDE.h"
#include "NodeDiagnostics.h"
#include "DiagReducer.h"
#include "Discretization.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::NodeDiagnostics;

NodeDiagnostics::NodeDiagnostics( const Discretization& d )
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization proxy to read from
// *****************************************************************************
{
  // Store the local IDs of those mesh nodes to which we contribute but do not
  // own, i.e., slave nodes. Ownership here is defined by having a lower chare
  // ID than any other chare that also contributes to the node.
  for (const auto& c : d.Msum())        // for all chares that neighbor our mesh
    if (d.thisIndex > c.first)          // if our chare ID is larger than theirs
      for (auto i : c.second)           // store local ID in set
        m_slave.insert( tk::cref_find( d.Lid(), i ) );
}

void
NodeDiagnostics::registerReducers()
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
NodeDiagnostics::compute( Discretization& d, const tk::Fields& u )
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization proxy to read from
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
    for (const auto& eq : g_cgpde)
      eq.initialize( d.Coord(), a, d.T()+d.Dt() );

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
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop(), 0.0 ) );

    // Put in norms sweeping our mesh chunk
    for (std::size_t i=0; i<u.nunk(); ++i)
      if (m_slave.find(i) == end(m_slave)) {    // ignore non-owned nodes
        // Compute sum for L2 norm of the numerical solution
        for (std::size_t c=0; c<u.nprop(); ++c)
          diag[L2SOL][c] += u(i,c,0) * u(i,c,0) * d.Vol()[i];
        // Compute sum for L2 norm of the numerical-analytic solution
        for (std::size_t c=0; c<u.nprop(); ++c)
          diag[L2ERR][c] +=
            (u(i,c,0)-a(i,c,0)) * (u(i,c,0)-a(i,c,0)) * d.Vol()[i];
        // Compute max for Linf norm of the numerical-analytic solution
        for (std::size_t c=0; c<u.nprop(); ++c) {
          auto err = std::abs( u(i,c,0) - a(i,c,0) );
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
