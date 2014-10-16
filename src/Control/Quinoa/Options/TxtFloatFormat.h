//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/TxtFloatFormat.h
  \author    J. Bakosi
  \date      Wed 01 Oct 2014 08:09:49 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Text floating-point output options and associations
  \details   Text floating-point output options and associations
*/
//******************************************************************************
#ifndef QuinoaTxtFloatFormatOptions_h
#define QuinoaTxtFloatFormatOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Txt floating-point format types
enum class TxtFloatFormatType : uint8_t { DEFAULT=0,
                                          FIXED,
                                          SCIENTIFIC };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, TxtFloatFormatType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class TxtFloatFormat : public tk::Toggle< TxtFloatFormatType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit TxtFloatFormat() :
      Toggle< TxtFloatFormatType >( "Text floating-point format",
      //! Enums -> names
      { { TxtFloatFormatType::DEFAULT, kw::txt_float_default().name() },
        { TxtFloatFormatType::FIXED, kw::txt_float_fixed().name() },
        { TxtFloatFormatType::SCIENTIFIC, kw::txt_float_scientific().name() } },
      //! keywords -> Enums
      { { kw::txt_float_default().string(), TxtFloatFormatType::DEFAULT },
        { kw::txt_float_fixed().string(), TxtFloatFormatType::FIXED },
        { kw::txt_float_scientific().string(), TxtFloatFormatType::SCIENTIFIC } }
      ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaTxtFloatFormatOptions_h
