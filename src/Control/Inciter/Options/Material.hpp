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
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Material types
enum class MaterialType : uint8_t { STIFFENEDGAS
                                  , JWL
                                  };

//! Pack/Unpack MaterialType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, MaterialType& e ) { PUP::pup( p, e ); }

//! \brief Material options: outsource to base templated on enum type
class Material : public tk::Toggle< MaterialType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::stiffenedgas
                                  , kw::jwl
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Material() :
      tk::Toggle< MaterialType >(
        //! Group, i.e., options, name
        kw::flux::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { MaterialType::STIFFENEDGAS, kw::stiffenedgas::name() }
        , { MaterialType::JWL, kw::jwl::name() }
        },
        //! keywords -> Enums
        { { kw::stiffenedgas::string(), MaterialType::STIFFENEDGAS }
        , { kw::jwl::string(), MaterialType::JWL }
        } )
    {}

};

} // ctr::
} // inciter::

#endif // MaterialOptions_h
