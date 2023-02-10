// *****************************************************************************
/*!
  \file      src/Control/Options/FieldFile.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field output file type options
  \details   Field output file type options
*/
// *****************************************************************************
#ifndef FieldFileOptions_h
#define FieldFileOptions_h

#include <brigand/sequences/list.hpp>

#include "QuinoaBuildConfig.hpp"
#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! Field output file types
enum class FieldFileType : uint8_t { EXODUSII
                                   };

//! \brief Pack/Unpack FieldFileType: forward overload to generic enum class
//!   packer
inline void operator|( PUP::er& p, FieldFileType& e ) { PUP::pup( p, e ); }

//! \brief FieldFileType options: outsource searches to base templated on enum
//!   type
class FieldFile : public tk::Toggle< FieldFileType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::exodusii
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit FieldFile() :
      tk::Toggle< FieldFileType >(
        //! Group, i.e., options, name 
        "Field output file type",
        //! Enums -> names
        { { FieldFileType::EXODUSII, kw::exodusii::name() },
        },
        //! keywords -> Enums
        { { kw::exodusii::string(), FieldFileType::EXODUSII },
        } ) {}
};

} // ctr::
} // tk:::

#endif // FieldFileOptions_h
