//******************************************************************************
/*!
  \file      src/Control/Options/PDFFile.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:08:39 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     PDF output file type options and associations
  \details   PDF output file type options and associations
*/
//******************************************************************************
#ifndef PDFFileOptions_h
#define PDFFileOptions_h

#include <map>

#include <Toggle.h>
#include <Keywords.h>

namespace tk {
namespace ctr {

//! PDF output file types
enum class PDFFileType : uint8_t { TXT=0,
                                   GMSHTXT,
                                   GMSHBIN,
                                   EXODUSII };

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
        { PDFFileType::GMSHTXT, kw::gmshtxt().name() },
        { PDFFileType::GMSHBIN, kw::gmshbin().name() },
        { PDFFileType::EXODUSII, kw::exodusii().name() } },
      //! keywords -> Enums
      { { kw::txt().string(), PDFFileType::TXT },
        { kw::gmshtxt().string(), PDFFileType::GMSHTXT },
        { kw::gmshbin().string(), PDFFileType::GMSHBIN },
        { kw::exodusii().string(), PDFFileType::EXODUSII } } ) {}
};

} // ctr::
} // tk:::

#endif // PDFFileOptions_h
