// *****************************************************************************
/*!
  \file      src/PDE/EoS/GetMatProp.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state class
  \details   This file defines functions for equations of state for the
    compressible flow equations.
*/
// *****************************************************************************
#ifndef GetMatProp_h
#define GetMatProp_h

#include <cmath>
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = tk::ncomp_t;

//! Get a property for a material
//! \tparam Prop Tag of property required
//! \param[in] imat Material-id who's property is required. Default is 0, so
//!   that for the single-material system, this argument can be left unspecified
//!   by the calling code
//! \return Material property Prop
//! \note This function returns a zero if the vector for the property required
//!   is empty. This will happen if the user has not specified that property
//!   in the control file, hence the inputdeck has not allocated that property
//!   vector.
template< class Prop >
tk::real
getmatprop( std::size_t imat=0 ) {
  const auto& matprop = g_inputdeck.get< tag::material >();
  const auto& map = g_inputdeck.get< tag::matidxmap >();
  auto meos = map.template get< tag::eosidx >()[ imat ];
  auto midx = map.template get< tag::matidx >()[ imat ];
  auto pvec = matprop[ meos ].template get< Prop >();

  tk::real mp;
  if (!pvec.empty())
    mp = pvec[ midx ];
  else
    mp = 0.0;
  return mp;
}

} //inciter::

#endif // GetMatProp_h
