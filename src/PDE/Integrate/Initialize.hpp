// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Initialize.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for initialization of system of PDEs in DG methods
  \details   This file contains functionality for setting initial conditions
     and evaluating known solutions used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Initialize_h
#define Initialize_h

#include "Basis.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

namespace tk {

//! Initalize a PDE system for DG by projecting the exact solution
//! in the DG solution space
void
initialize( ncomp_t system,
            ncomp_t ncomp,
            ncomp_t offset,
            const std::vector< inciter::EoS_Base* >& mat_blk,
            const Fields& L,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const InitializeFn& solution,
            Fields& unk,
            real t,
            const std::size_t nielem );

//! Update the rhs by adding the initial analytical solution term
void
update_rhs( ncomp_t ncomp,
            const std::size_t ndof,
            const tk::real wt,
            const std::vector< tk::real >& B,
            const std::vector< tk::real >& s,
            std::vector< tk::real >& R );


//! Compute the initial conditions
void
eval_init( ncomp_t ncomp,
           ncomp_t offset,
           const std::size_t ndof,
           const std::size_t rdof,
           const std::size_t e,
           const std::vector< tk::real >& R,
           const Fields& L,
           Fields& unk );

template< class eq >
void
BoxElems( std::size_t system,
  const tk::Fields& geoElem,
  std::size_t nielem,
  std::vector< std::unordered_set< std::size_t > >& inbox )
// *****************************************************************************
//! Determine if elements lie inside user defined IC boxes
//! \tparam eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] geoElem Element geometry array
//! \param[in] nielem Number of internal elements
//! \param[in,out] inbox List of nodes at which box user ICs are set for
//!    each IC box
// *****************************************************************************
{
  // Detect if user has configured a IC boxes
  const auto& icbox = inciter::g_inputdeck.get<tag::param, eq, tag::ic,
    tag::box>();
  if (icbox.size() > system) {
    std::size_t bcnt = 0;
    for (const auto& b : icbox[system]) {   // for all boxes for this eq
     inbox.emplace_back();
      std::vector< tk::real > box
        { b.template get< tag::xmin >(), b.template get< tag::xmax >(),
          b.template get< tag::ymin >(), b.template get< tag::ymax >(),
          b.template get< tag::zmin >(), b.template get< tag::zmax >() };

      const auto eps = std::numeric_limits< tk::real >::epsilon();
      // Determine which elements lie in the IC box
      for (ncomp_t e=0; e<nielem; ++e) {
        auto x = geoElem(e,1,0);
        auto y = geoElem(e,2,0);
        auto z = geoElem(e,3,0);
        if ( std::any_of( begin(box), end(box),
                          [=](auto p){ return abs(p) > eps; } ) &&
             x>box[0] && x<box[1] &&
             y>box[2] && y<box[3] &&
             z>box[4] && z<box[5] )
        {
          inbox[bcnt].insert( e );
        }
      }
      ++bcnt;
    }
  }
}

} // tk::

#endif // Initialize_h
