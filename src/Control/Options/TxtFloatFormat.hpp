// *****************************************************************************
/*!
  \file      src/Control/Options/TxtFloatFormat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Text floating-point output options
  \details   Text floating-point output options
*/
// *****************************************************************************
#ifndef TxtFloatFormatOptions_h
#define TxtFloatFormatOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! Txt floating-point format types
enum class TxtFloatFormatType : uint8_t { DEFAULT=0,
                                          FIXED,
                                          SCIENTIFIC };

//! \brief Pack/Unpack TxtFloatFormatType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, TxtFloatFormatType& e ) { PUP::pup( p, e ); }

//! \brief TxtFloatFormat options: outsource searches to base templated on enum
//!   type
class TxtFloatFormat : public tk::Toggle< TxtFloatFormatType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::txt_float_default
                                  , kw::txt_float_fixed
                                  , kw::txt_float_scientific
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit TxtFloatFormat() :
      tk::Toggle< TxtFloatFormatType >(
        //! Group, i.e., options, name
        "floating-point format",
        //! Enums -> names
        { { TxtFloatFormatType::DEFAULT, kw::txt_float_default::name() },
          { TxtFloatFormatType::FIXED, kw::txt_float_fixed::name() },
          { TxtFloatFormatType::SCIENTIFIC, kw::txt_float_scientific::name() } },
        //! keywords -> Enums
        { { kw::txt_float_default::string(), TxtFloatFormatType::DEFAULT },
          { kw::txt_float_fixed::string(), TxtFloatFormatType::FIXED },
        { kw::txt_float_scientific::string(), TxtFloatFormatType::SCIENTIFIC } }
      ) {}
};

} // ctr::
} // tk:::

#endif // TxtFloatFormatOptions_h
