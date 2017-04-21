// *****************************************************************************
/*!
  \file      src/Control/Options/TxtFloatFormat.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Text floating-point output options
  \details   Text floating-point output options
*/
// *****************************************************************************
#ifndef TxtFloatFormatOptions_h
#define TxtFloatFormatOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! Txt floating-point format types
//! \author J. Bakosi
enum class TxtFloatFormatType : uint8_t { DEFAULT=0,
                                          FIXED,
                                          SCIENTIFIC };

//! \brief Pack/Unpack TxtFloatFormatType: forward overload to generic enum
//!   class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, TxtFloatFormatType& e ) { PUP::pup( p, e ); }

//! \brief TxtFloatFormat options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class TxtFloatFormat : public tk::Toggle< TxtFloatFormatType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::txt_float_default
                                       , kw::txt_float_fixed
                                       , kw::txt_float_scientific
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
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
