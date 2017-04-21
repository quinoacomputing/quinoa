// *****************************************************************************
/*!
  \file      src/Control/Options/RNGSSESeqLen.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     RNGSSE sequence length options
  \details   RNGSSE sequence length options
*/
// *****************************************************************************
#ifndef RNGSSESeqLenOptions_h
#define RNGSSESeqLenOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! RNGSSE sequence length options
//! \author J. Bakosi
enum class RNGSSESeqLenType : uint8_t { SHORT,
                                        MEDIUM,
                                        LONG };

//! \brief Pack/Unpack RNGSSESeqLenType: forward overload to generic enum class
//!   packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, RNGSSESeqLenType& e ) { PUP::pup( p, e ); }

//! \brief RNGSSESeqLen options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class RNGSSESeqLen : public tk::Toggle< RNGSSESeqLenType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::seq_short
                                       , kw::seq_med
                                       , kw::seq_long
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit RNGSSESeqLen() :
      tk::Toggle< RNGSSESeqLenType >(
        //! Group, i.e., options, name
        "sequence length",
        //! Enums -> names
        { { RNGSSESeqLenType::SHORT, kw::seq_short::name() },
          { RNGSSESeqLenType::MEDIUM, kw::seq_med::name() },
          { RNGSSESeqLenType::LONG, kw::seq_long::name() } },
        //! keywords -> Enums
        { { kw::seq_short::string(), RNGSSESeqLenType::SHORT },
          { kw::seq_med::string(), RNGSSESeqLenType::MEDIUM },
          { kw::seq_long::string(), RNGSSESeqLenType::LONG } } ) {}
};

} // ctr::
} // tk::

#endif // RNGSSESeqLenOptions_h
