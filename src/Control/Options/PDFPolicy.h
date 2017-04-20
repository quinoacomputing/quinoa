// *****************************************************************************
/*!
  \file      src/Control/Options/PDFPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     PDF output file policy options
  \details   PDF output file policy options
*/
// *****************************************************************************
#ifndef PDFPolicyOptions_h
#define PDFPolicyOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! PDF output file policy types
//! \author J. Bakosi
enum class PDFPolicyType : uint8_t { OVERWRITE=0,
                                     MULTIPLE,
                                     EVOLUTION };

//! \brief Pack/Unpack PDFPolicyType: forward overload to generic enum class
//!   packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, PDFPolicyType& e ) { PUP::pup( p, e ); }

//! \brief PDFPolicy ptions: outsource searches to base templated on enum type
//! \author J. Bakosi
class PDFPolicy : public tk::Toggle< PDFPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::overwrite
                                       , kw::multiple
                                       , kw::evolution
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
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
