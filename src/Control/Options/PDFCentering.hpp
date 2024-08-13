// *****************************************************************************
/*!
  \file      src/Control/Options/PDFCentering.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     PDF output file centering type options
  \details   PDF output file centering type options
*/
// *****************************************************************************
#ifndef PDFCenteringOptions_h
#define PDFCenteringOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! PDF output file types
enum class PDFCenteringType : uint8_t { ELEM=0,
                                        NODE };

//! \brief Pack/Unpack PDFCenteringType: forward overload to generic enum class
//!   packer
inline void operator|( PUP::er& p, PDFCenteringType& e ) { PUP::pup( p, e ); }

//! \brief PDFCentering options: outsource searches to base templated on enum
//!   type
class PDFCentering : public tk::Toggle< PDFCenteringType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PDFCentering() :
      tk::Toggle< PDFCenteringType >(
        //! Group, i.e., options, name
        "PDF output file centering",
        //! Enums -> names
        { { PDFCenteringType::ELEM, "elem" },
          { PDFCenteringType::NODE, "node" } },
        //! keywords -> Enums
        { { "elem", PDFCenteringType::ELEM },
          { "node", PDFCenteringType::NODE } } ) {}
};

} // ctr::
} // tk:::

#endif // PDFCenteringOptions_h
