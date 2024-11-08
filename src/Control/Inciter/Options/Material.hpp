// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Material.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Material types for inciter
  \details   Material types for inciter.
*/
// *****************************************************************************
#ifndef MaterialOptions_h
#define MaterialOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Material types
enum class MaterialType : uint8_t { STIFFENEDGAS
                                  , JWL
                                  , SMALLSHEARSOLID
                                  , GODUNOVROMENSKIALUMINUM
                                  , THERMALLYPERFECTGAS
                                  };

//! Pack/Unpack MaterialType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, MaterialType& e ) { PUP::pup( p, e ); }

//! \brief Material options: outsource to base templated on enum type
class Material : public tk::Toggle< MaterialType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Material() :
      tk::Toggle< MaterialType >(
        //! Group, i.e., options, name
        "Material EOS",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { MaterialType::STIFFENEDGAS, "stiffenedgas" }
        , { MaterialType::JWL, "jwl" }
        , { MaterialType::SMALLSHEARSOLID, "smallshearsolid" }
        , { MaterialType::GODUNOVROMENSKIALUMINUM, "godunovromenski_aluminum" }
        , { MaterialType::THERMALLYPERFECTGAS, "thermallyperfectgas" }
        },
        //! keywords -> Enums
        { { "stiffenedgas", MaterialType::STIFFENEDGAS }
        , { "jwl", MaterialType::JWL }
        , { "smallshearsolid", MaterialType::SMALLSHEARSOLID }
        , { "godunovromenski_aluminum", MaterialType::GODUNOVROMENSKIALUMINUM }
        , { "thermallyperfectgas", MaterialType::THERMALLYPERFECTGAS }
        } )
    {}

};

} // ctr::
} // inciter::

#endif // MaterialOptions_h
