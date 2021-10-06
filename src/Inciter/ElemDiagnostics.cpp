// *****************************************************************************
/*!
  \file      src/Inciter/ElemDiagnostics.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ElemDiagnostics class for collecting element diagnostics
  \details   ElemDiagnostics class for collecting element diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <cmath>

#include "DGPDE.hpp"
#include "ElemDiagnostics.hpp"
#include "DiagReducer.hpp"
#include "Discretization.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

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
ElemDiagnostics::compute( Discretization& d,
                          const std::size_t nchGhost,
                          const tk::Fields& geoElem,
                          const std::vector< std::size_t >& ndofel,
                          const tk::Fields& u ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
//! \param[in] d Discretization base class to read from
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] ndofel Vector of local number of degrees of freedom
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
  auto diagfreq = g_inputdeck.get< tag::output, tag::iter, tag::diag >();

  if ( !((d.It()+1) % diagfreq) ) {  // if remainder, don't compute diagnostics

    // Query number of degrees of freedom from user's setting
    const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

    // Diagnostics vector (of vectors) during aggregation. See
    // Inciter/Diagnostics.h.
    std::vector< std::vector< tk::real > >
      diag( NUMDIAG, std::vector< tk::real >( u.nprop()/rdof, 0.0 ) );

    // Compute diagnostics for DG
    compute_diag(d, rdof, nchGhost, geoElem, ndofel, u, diag);

    // Append diagnostics vector with metadata on the current time step
    // ITER: Current iteration count (only the first entry is used)
    // TIME: Current physical time (only the first entry is used)
    // DT: Current physical time step size (only the first entry is used)
    diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
    diag[TIME][0] = d.T() + d.Dt();
    diag[DT][0] = d.Dt();

    // Contribute to diagnostics
    auto stream = serialize( d.MeshId(), diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;        // diagnostics have been computed

  }

  return false;         // diagnostics have not been computed
}

void
ElemDiagnostics::compute_diag( const Discretization& d,
                               const std::size_t rdof,
                               const std::size_t nchGhost,
                               const tk::Fields& geoElem,
                               const std::vector< std::size_t >& ndofel,
                               const tk::Fields& u,
                               std::vector< std::vector< tk::real > >& diag )
const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc. for DG
//! \param[in] d Discretization base class to read from
//! \param[in] rdof Number of reconstructed degrees of freedom
//! \param[in] nchGhost Number of chare boundary ghost elements
//! \param[in] geoElem Element geometry
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] u Current solution vector
//! \param[in,out] diag Diagnostics vector
// *****************************************************************************
{
  const auto& inpoel = d.Inpoel();
  const auto& coord = d.Coord();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<u.nunk()-nchGhost; ++e)
  {
    // Number of quadrature points for volume integration
    auto ng = tk::NGdiag(ndofel[e]);

    // arrays for quadrature points
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = tk::eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B = tk::eval_basis( ndofel[e], coordgp[0][igp], coordgp[1][igp],
                               coordgp[2][igp]);

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      std::vector< tk::real > s;

      for (const auto& eq : g_dgpde)
        // cppcheck-suppress useStlAlgorithm
        s = eq.solution( gp[0], gp[1], gp[2], d.T()+d.Dt() );

      for (std::size_t c=0; c<u.nprop()/rdof; ++c)
      {
        auto mark = c*rdof;
        auto ugp = u(e, mark, 0);

        if(ndofel[e] > 1)
        {
          ugp +=  u(e, mark+1, 0) * B[1]
                + u(e, mark+2, 0) * B[2]
                + u(e, mark+3, 0) * B[3];

          if(ndofel[e] > 4)
          {
            ugp +=  u(e, mark+4, 0) * B[4]
                  + u(e, mark+5, 0) * B[5]
                  + u(e, mark+6, 0) * B[6]
                  + u(e, mark+7, 0) * B[7]
                  + u(e, mark+8, 0) * B[8]
                  + u(e, mark+9, 0) * B[9];
          }
        }

        // Compute sum for L2 norm of the numerical solution
        diag[L2SOL][c] += wt * ugp * ugp;

        // Compute sum for L2 norm of the numerical-analytic solution
        diag[L2ERR][c] += wt * (ugp-s[c]) * (ugp-s[c]);

        // Compute max for Linf norm of the numerical-analytic solution
        auto err = std::abs( ugp - s[c] );
        if (err > diag[LINFERR][c]) diag[LINFERR][c] = err;

        // Compute sum of the total energy over the entire domain (only the
        // first entry is used)
        if (c == u.nprop()/rdof-1) diag[TOTALSOL][0] += wt * ugp;
      }
    }
  }
}
