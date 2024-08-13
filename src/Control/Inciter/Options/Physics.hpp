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
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Physics types
enum class PhysicsType : uint8_t { ADVECTION
                                 , ADVDIFF
                                 , EULER
                                 , ENERGYPILL
                                 };

//! Pack/Unpack PhysicsType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PhysicsType& e ) { PUP::pup( p, e ); }

//! \brief Physics options: outsource to base templated on enum type
class Physics : public tk::Toggle< PhysicsType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Physics() :
      tk::Toggle< PhysicsType >(
        //! Group, i.e., options, name
        "Physics",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { PhysicsType::ADVECTION, "advection" }
        , { PhysicsType::ADVDIFF, "advdiff" }
        , { PhysicsType::EULER, "euler" }
        , { PhysicsType::ENERGYPILL, "energy_pill" }
        },
        //! keywords -> Enums
        { { "advection", PhysicsType::ADVECTION }
        , { "advdiff", PhysicsType::ADVDIFF }
        , { "euler", PhysicsType::EULER }
        , { "energy_pill", PhysicsType::ENERGYPILL }
        } )
    {}
};

} // ctr::
} // inciter::

#endif // PhysicsOptions_h
