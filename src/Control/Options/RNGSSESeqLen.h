//******************************************************************************
/*!
  \file      src/Control/Options/RNGSSESeqLen.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:08:58 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGSSE sequence length options
  \details   RNGSSE sequence length options
*/
//******************************************************************************
#ifndef RNGSSESeqLenOptions_h
#define RNGSSESeqLenOptions_h

#include <map>

#include <Toggle.h>
#include <Keywords.h>

namespace tk {
namespace ctr {

//! RNGSSE's sequence length options
enum class RNGSSESeqLenType : uint8_t { SHORT,
                                        MEDIUM,
                                        LONG };

//! Pack/Unpack: delegate to tk::
inline void operator|( PUP::er& p, RNGSSESeqLenType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class RNGSSESeqLen : public tk::Toggle< RNGSSESeqLenType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNGSSESeqLen() :
      Toggle< RNGSSESeqLenType >( "sequence length",
        //! Enums -> names
        { { RNGSSESeqLenType::SHORT, kw::seq_short().name() },
          { RNGSSESeqLenType::MEDIUM, kw::seq_med().name() },
          { RNGSSESeqLenType::LONG, kw::seq_long().name() } },
        //! keywords -> Enums
        { { kw::seq_short().string(), RNGSSESeqLenType::SHORT },
          { kw::seq_med().string(), RNGSSESeqLenType::MEDIUM },
          { kw::seq_long().string(), RNGSSESeqLenType::LONG } } ) {}
};

} // ctr::
} // tk::

#endif // RNGSSESeqLenOptions_h
