// *****************************************************************************
/*!
  \file      src/Control/Options/PDFPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     PDF output file policy options
  \details   PDF output file policy options
*/
// *****************************************************************************
#ifndef PDFPolicyOptions_h
#define PDFPolicyOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! PDF output file policy types
enum class PDFPolicyType : uint8_t { OVERWRITE=0,
                                     MULTIPLE,
                                     EVOLUTION };

//! \brief Pack/Unpack PDFPolicyType: forward overload to generic enum class
//!   packer
inline void operator|( PUP::er& p, PDFPolicyType& e ) { PUP::pup( p, e ); }

//! \brief PDFPolicy ptions: outsource searches to base templated on enum type
class PDFPolicy : public tk::Toggle< PDFPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::overwrite
                                  , kw::multiple
                                  , kw::evolution
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PDFPolicy() :
      tk::Toggle< PDFPolicyType >(
        //! Group, i.e., options, name
        "PDF output file policy",
        //! Enums -> names
        { { PDFPolicyType::OVERWRITE, kw::overwrite::name() },
          { PDFPolicyType::MULTIPLE, kw::multiple::name() },
          { PDFPolicyType::EVOLUTION, kw::evolution::name() } },
        //! keywords -> Enums
        { { kw::overwrite::string(), PDFPolicyType::OVERWRITE },
          { kw::multiple::string(), PDFPolicyType::MULTIPLE },
          { kw::evolution::string(), PDFPolicyType::EVOLUTION } } ) {}
};

} // ctr::
} // tk:::

#endif // PDFPolicyOptions_h
