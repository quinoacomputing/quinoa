//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/PDFCentering.h
  \author    J. Bakosi
  \date      Wed 17 Sep 2014 10:31:34 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     PDF output file centering type options and associations
  \details   PDF output file centering type options and associations
*/
//******************************************************************************
#ifndef QuinoaPDFCenteringOptions_h
#define QuinoaPDFCenteringOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! PDF output file types
enum class PDFCenteringType : uint8_t { ELEM=0,
                                        NODE };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDFCenteringType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class PDFCentering : public tk::Toggle< PDFCenteringType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDFCentering() :
      Toggle< PDFCenteringType >( "PDF output file centering",
      //! Enums -> names
      { { PDFCenteringType::ELEM, kw::elem().name() },
        { PDFCenteringType::NODE, kw::node().name() } },
      //! keywords -> Enums
      { { kw::elem().string(), PDFCenteringType::ELEM },
        { kw::node().string(), PDFCenteringType::NODE } } ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaPDFCenteringOptions_h
