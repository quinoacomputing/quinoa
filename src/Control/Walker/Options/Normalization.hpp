// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/Normalization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Normalization variants (used for MixDirichlet)
  \details   Normalization variants (used for MixDirichlet).
*/
// *****************************************************************************
#ifndef NormalizationOptions_h
#define NormalizationOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace walker {
namespace ctr {

//! Normalization types
enum class NormalizationType : uint8_t { LIGHT=0
                                       , HEAVY };

//! \brief Pack/Unpack NormalizationType: forward overload to generic enum
//!    class packer
inline void operator|( PUP::er& p, NormalizationType& e )
{ PUP::pup( p, e ); }

//! Normalization variants: outsource to base templated on enum type
class Normalization : public tk::Toggle< NormalizationType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::light
                                  , kw::heavy
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Normalization() :
      tk::Toggle< NormalizationType >(
        //! Group, i.e., options, name
        "Normalization",
        //! Enums -> names
        { { NormalizationType::LIGHT, kw::light::name() },
          { NormalizationType::HEAVY, kw::heavy::name() } },
        //! keywords -> Enums
        { { kw::light::string(), NormalizationType::LIGHT },
          { kw::heavy::string(), NormalizationType::HEAVY } } )
    {}
};

} // ctr::
} // walker::

#endif // NormalizationOptions_h
