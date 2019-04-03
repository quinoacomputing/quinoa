// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/BC.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary condition options for inciter
  \details   Boundary condition options for inciter
*/
// *****************************************************************************
#ifndef InciterBCOptions_h
#define InciterBCOptions_h

#include <brigand/sequences/list.hpp>

#include "TaggedTuple.hpp"
#include "Toggle.hpp"
#include "Keywords.hpp"

namespace inciter {
namespace ctr {

//! Boundary condition types
enum class BCType : uint8_t { SYM,
                              INLET,
                              OUTLET,
                              EXTRAPOLATE };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, BCType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class BC : public tk::Toggle< BCType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::bc_sym
                                  , kw::bc_inlet
                                  , kw::bc_outlet
                                  , kw::bc_extrapolate
                                  >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit BC() :
      tk::Toggle< BCType >( "Boundary condition",
        //! Enums -> names
        { { BCType::SYM, kw::bc_sym::name() },
          { BCType::INLET, kw::bc_inlet::name() },
          { BCType::OUTLET, kw::bc_outlet::name() },
          { BCType::EXTRAPOLATE, kw::bc_extrapolate::name() } },
        //! keywords -> Enums
        { { kw::bc_sym::string(), BCType::SYM },
          { kw::bc_inlet::string(), BCType::INLET },
          { kw::bc_outlet::string(), BCType::OUTLET },
          { kw::bc_extrapolate::string(), BCType::EXTRAPOLATE } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterBCOptions_h
