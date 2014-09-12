//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/PDFFile.h
  \author    J. Bakosi
  \date      Wed 10 Sep 2014 05:28:22 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     PDFFile model options and associations
  \details   PDFFile model options and associations
*/
//******************************************************************************
#ifndef QuinoaPDFFileOptions_h
#define QuinoaPDFFileOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! PDFFile model types
enum class PDFFileType : uint8_t { OVERWRITE=0,
                                   MULTIPLE };

//! Class with base templated on the above enum class with associations
class PDFFile : public tk::Toggle< PDFFileType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDFFile() :
      Toggle< PDFFileType >( "PDF file type",
      //! Enums -> names
      { { PDFFileType::OVERWRITE, kw::overwrite().name() },
        { PDFFileType::MULTIPLE, kw::multiple().name() } },
      //! keywords -> Enums
      { { kw::overwrite().string(), PDFFileType::OVERWRITE },
        { kw::multiple().string(), PDFFileType::MULTIPLE } } ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaPDFFileOptions_h
