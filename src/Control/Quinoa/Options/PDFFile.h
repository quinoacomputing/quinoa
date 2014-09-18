//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/PDFFile.h
  \author    J. Bakosi
  \date      Tue 16 Sep 2014 08:14:53 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     PDF output file type options and associations
  \details   PDF output file type options and associations
*/
//******************************************************************************
#ifndef QuinoaPDFFileOptions_h
#define QuinoaPDFFileOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! PDF output file types
enum class PDFFileType : uint8_t { TXT=0,
                                   GMSH };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDFFileType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class PDFFile : public tk::Toggle< PDFFileType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDFFile() :
      Toggle< PDFFileType >( "PDF output file type",
      //! Enums -> names
      { { PDFFileType::TXT, kw::txt().name() },
        { PDFFileType::GMSH, kw::gmsh().name() } },
      //! keywords -> Enums
      { { kw::txt().string(), PDFFileType::TXT },
        { kw::gmsh().string(), PDFFileType::GMSH } } ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaPDFFileOptions_h
