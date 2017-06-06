// *****************************************************************************
/*!
  \file      src/Control/Options/Error.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Error type options
  \details   Error type options
*/
// *****************************************************************************
#ifndef ErrorOptions_h
#define ErrorOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! Error types
//! \author J. Bakosi
enum class ErrorType : uint8_t { L2=0,
                                 LINF };

//! \brief Pack/Unpack ErrorType: forward overload to generic enum
//!   class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, ErrorType& e ) { PUP::pup( p, e ); }

//! \brief Error options: outsource searches to base templated on enum type
//! \author J. Bakosi
class Error : public tk::Toggle< ErrorType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::l2
                                       , kw::linf
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Error() :
      tk::Toggle< ErrorType >(
        //! Group, i.e., options, name
        "error",
        //! Enums -> names
        { { ErrorType::L2, kw::l2::name() },
          { ErrorType::LINF, kw::linf::name() } },
        //! keywords -> Enums
        { { kw::l2::string(), ErrorType::L2 },
          { kw::linf::string(), ErrorType::LINF } }
      ) {}
};

} // ctr::
} // tk:::

#endif // ErrorOptions_h
