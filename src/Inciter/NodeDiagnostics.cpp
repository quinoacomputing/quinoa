// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     NodeDiagnostics class for collecting nodal diagnostics
  \details   NodeDiagnostics class for collecting nodal diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************

#include "CGPDE.hpp"
#include "NodeDiagnostics.hpp"
#include "DiagReducer.hpp"
#include "Discretization.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Refiner.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::NodeDiagnostics;

void
NodeDiagnostics::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details This routine is supposed to be called from a Charm++ initnode
//!   routine. Since the runtime system executes initnode routines exactly once
//!   on every logical node early on in the Charm++ init sequence, they must be
//!   static as they are called without an object. See also: Section
//!   "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  DiagMerger = CkReduction::addReducer( mergeDiag );
}

bool
NodeDiagnostics::compute(
  Discretization& d,
  const tk::Fields& u,
  const tk::Fields& un,
  const std::unordered_map< int,
          std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
  const std::unordered_set< std::size_t >& symbcnodes,
  const std::unordered_set< std::size_t >& farfieldbcnodes ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization proxy to read from
//! \param[in] u Current solution vector
//! \param[in] un Previous solution vector
//! \param[in] bnorm Face normals in boundary points, key local node id,
//!   first 3 reals of value: unit normal, outer key: side set id
//! \param[in] symbcnodes Unique set of node ids at which to set symmetry BCs
//! \param[in] farfieldbcnodes Unique set of node ids at which to set farfield
//!   BCs
//! \return True if diagnostics have been computed
//! \details Diagnostics are defined as some norm, e.g., L2 norm, of a quantity,
//!   computed in mesh nodes, A, as ||A||_2 = sqrt[ sum_i(A_i)^2 V_i ],
//!   where the sum is taken over all mesh nodes and V_i is the nodal volume.
//!   We send multiple sets of quantities to the host for aggregation across
//!   the whole mesh. The final aggregated solution will end up in
//!   Transporter::diagnostics(). Aggregation of the partially computed
//!   diagnostics is done via potentially different policies for each field.
//! \see inciter::mergeDiag(), src/Inciter/Diagnostics.hpp
// *****************************************************************************
{
  // Optionally collect diagnostics and send for aggregation across all workers

  // Query after how many time steps user wants to dump diagnostics
  auto diagfreq = g_inputdeck.get< tag::diagnostics, tag::interval >();

  if ( !((d.It()+1) % diagfreq) ) {     // if remainder, don't dump

    // Diagnostics vector (of vectors) during aggregation. See
    // Inciter/Diagnostics.h.
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop(), 0.0 ) );

    const auto& coord = d.Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    const auto& v = d.V();  // nodal volumes without contributions from others

    // Evaluate analytic solution (if exist, if not, IC)
    auto an = u;
    auto mv = d.meshvel();
    for (std::size_t i=0; i<an.nunk(); ++i) {
      // Query analytic solution for all components of all PDEs integrated
      std::vector< tk::real > a;
      auto s = g_cgpde[d.MeshId()].solution( x[i], y[i], z[i], d.T()+d.Dt() );
      std::move( begin(s), end(s), std::back_inserter(a) );
      Assert( a.size() == u.nprop(), "Size mismatch" );
      for (std::size_t c=0; c<an.nprop(); ++c) an(i,c) = a[c];
    }
    // Apply symmetry BCs on analytic solution (if exist, if not, IC)
    g_cgpde[d.MeshId()].symbc( an, mv, coord, bnorm, symbcnodes );
    // Apply farfield BCs on analytic solution (if exist, if not, IC)
    g_cgpde[d.MeshId()].farfieldbc( an, coord, bnorm, farfieldbcnodes );

    // Put in norms sweeping our mesh chunk
    for (std::size_t i=0; i<u.nunk(); ++i) {
      // Compute sum for L2 norm of the numerical solution
      for (std::size_t c=0; c<u.nprop(); ++c)
        diag[L2SOL][c] += u(i,c) * u(i,c) * v[i];
      // Compute sum for L2 norm of the numerical-analytic solution
      for (std::size_t c=0; c<u.nprop(); ++c)
        diag[L2ERR][c] += (u(i,c)-an(i,c)) * (u(i,c)-an(i,c)) * v[i];
      // Compute sum for L2 norm of the residual
      for (std::size_t c=0; c<u.nprop(); ++c)
        diag[L2RES][c] += (u(i,c)-un(i,c)) * (u(i,c)-un(i,c)) * v[i];
      // Compute max for Linf norm of the numerical-analytic solution
      for (std::size_t c=0; c<u.nprop(); ++c) {
        auto err = std::abs( u(i,c) - an(i,c) );
        if (err > diag[LINFERR][c]) diag[LINFERR][c] = err;
      }
      // Compute sum of the total energy over the entire domain (only the first
      // entry is used)
      diag[TOTALSOL][0] += u(i,u.nprop()-1) * v[i];
    }

    // Append diagnostics vector with metadata on the current time step
    // ITER:: Current iteration count (only the first entry is used)
    // TIME: Current physical time (only the first entry is used)
    // DT: Current physical time step size (only the first entry is used)
    diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
    diag[TIME][0] = d.T() + d.Dt();
    diag[DT][0] = d.Dt();

    // Contribute to diagnostics
    auto stream = serialize( d.MeshId(), u.nprop(), diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;        // diagnostics have been computed
  }

  return false;         // diagnostics have not been computed
}
