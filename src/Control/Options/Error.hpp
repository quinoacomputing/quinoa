// *****************************************************************************
/*!
  \file      src/Control/Options/Error.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Error type options
  \details   Error type options
*/
// *****************************************************************************
#ifndef ErrorOptions_h
#define ErrorOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! Error types
enum class ErrorType : uint8_t { L2=0,
                                 LINF };

//! \brief Pack/Unpack ErrorType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, ErrorType& e ) { PUP::pup( p, e ); }

//! \brief Error options: outsource searches to base templated on enum type
class Error : public tk::Toggle< ErrorType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Error() :
      tk::Toggle< ErrorType >(
        //! Group, i.e., options, name
        "error",
        //! Enums -> names
        { { ErrorType::L2, "l2" },
          { ErrorType::LINF, "linf" } },
        //! keywords -> Enums
        { { "l2", ErrorType::L2 },
          { "linf", ErrorType::LINF } }
      ) {}
};

} // ctr::
} // tk:::

#endif // ErrorOptions_h
