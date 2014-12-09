//******************************************************************************
/*!
  \file      src/Control/Options/PDFPolicy.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:08:45 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     PDF output file policy options and associations
  \details   PDF output file policy options and associations
*/
//******************************************************************************
#ifndef PDFPolicyOptions_h
#define PDFPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Keywords.h>

namespace tk {
namespace ctr {

//! PDF output file policy types
enum class PDFPolicyType : uint8_t { OVERWRITE=0,
                                     MULTIPLE,
                                     EVOLUTION };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDFPolicyType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class PDFPolicy : public tk::Toggle< PDFPolicyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDFPolicy() :
      Toggle< PDFPolicyType >( "PDF output file policy",
      //! Enums -> names
      { { PDFPolicyType::OVERWRITE, kw::overwrite().name() },
        { PDFPolicyType::MULTIPLE, kw::multiple().name() },
        { PDFPolicyType::EVOLUTION, kw::evolution().name() } },
      //! keywords -> Enums
      { { kw::overwrite().string(), PDFPolicyType::OVERWRITE },
        { kw::multiple().string(), PDFPolicyType::MULTIPLE },
        { kw::evolution().string(), PDFPolicyType::EVOLUTION } } ) {}
};

} // ctr::
} // walker:::

#endif // PDFPolicyOptions_h
