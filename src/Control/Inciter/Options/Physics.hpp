// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Physics.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics options for inciter
  \details   Physics options for inciter
*/
// *****************************************************************************
#ifndef PhysicsOptions_h
#define PhysicsOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Physics types
enum class PhysicsType : uint8_t { ADVECTION,
                                   ADVDIFF,
                                   EULER,
                                   NAVIERSTOKES,
                                   ENERGYPILL };

//! Pack/Unpack PhysicsType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PhysicsType& e ) { PUP::pup( p, e ); }

//! \brief Physics options: outsource to base templated on enum type
class Physics : public tk::Toggle< PhysicsType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::advection
                                  , kw::advdiff
                                  , kw::euler
                                  , kw::navierstokes
                                  , kw::energy_pill
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Physics() :
      tk::Toggle< PhysicsType >(
        //! Group, i.e., options, name
        kw::physics::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { PhysicsType::ADVECTION, kw::advection::name() },
          { PhysicsType::ADVDIFF, kw::advdiff::name() },
          { PhysicsType::EULER, kw::euler::name() },
          { PhysicsType::NAVIERSTOKES, kw::navierstokes::name() },
          { PhysicsType::ENERGYPILL, kw::energy_pill::name() }
        },
        //! keywords -> Enums
        { { kw::advection::string(), PhysicsType::ADVECTION },
          { kw::advdiff::string(), PhysicsType::ADVDIFF },
          { kw::euler::string(), PhysicsType::EULER },
          { kw::navierstokes::string(), PhysicsType::NAVIERSTOKES },
          { kw::energy_pill::string(), PhysicsType::ENERGYPILL }
        } )
    {}
};

} // ctr::
} // inciter::

#endif // PhysicsOptions_h
