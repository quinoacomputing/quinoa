// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/MeshVelocity.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh velocity configuration options for inciter
  \details   Mesh velocity configuration options for inciter.
*/
// *****************************************************************************
#ifndef MeshVelocityOptions_h
#define MeshVelocityOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Mesh velocity configuration option types
enum class MeshVelocityType : uint8_t { SINE
                                      , FLUID
                                      , USER_DEFINED
                                      };

//! Pack/Unpack MeshVelocityType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, MeshVelocityType& e ) { PUP::pup( p, e ); }

//! \brief Mesh velocity options: outsource to base templated on enum type
class MeshVelocity : public tk::Toggle< MeshVelocityType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::sine
                                  , kw::fluid
                                  , kw::user_defined
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MeshVelocity() :
      tk::Toggle< MeshVelocityType >(
        //! Group, i.e., options, name
        kw::meshvelocity::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { MeshVelocityType::SINE, kw::sine::name() }
        , { MeshVelocityType::FLUID, kw::fluid::name() }
        , { MeshVelocityType::USER_DEFINED, kw::user_defined::name() }
        },
        //! keywords -> Enums
        { { kw::sine::string(), MeshVelocityType::SINE }
        , { kw::fluid::string(), MeshVelocityType::FLUID }
        , { kw::user_defined::string(), MeshVelocityType::USER_DEFINED }
        } )
    {}

};

} // ctr::
} // inciter::

#endif // MeshVelocityOptions_h
