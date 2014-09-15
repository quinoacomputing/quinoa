//******************************************************************************
/*!
  \file      src/Control/Options/RNGSSESeqLen.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 08:02:00 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNGSSE sequence length options
  \details   RNGSSE sequence length options
*/
//******************************************************************************
#ifndef RNGSSESeqLenOptions_h
#define RNGSSESeqLenOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace tk {
namespace ctr {

//! RNGSSE's sequence length options
enum class RNGSSESeqLenType : uint8_t { SHORT,
                                        MEDIUM,
                                        LONG };

//! Pack/Unpack: delegate to tk::
inline void operator|( PUP::er& p, RNGSSESeqLenType& e ) { tk::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class RNGSSESeqLen : public tk::Toggle< RNGSSESeqLenType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNGSSESeqLen() :
      Toggle< RNGSSESeqLenType >( "sequence length",
        //! Enums -> names
        { { RNGSSESeqLenType::SHORT, tk::kw::seq_short().name() },
          { RNGSSESeqLenType::MEDIUM, tk::kw::seq_med().name() },
          { RNGSSESeqLenType::LONG, tk::kw::seq_long().name() } },
        //! keywords -> Enums
        { { tk::kw::seq_short().string(), RNGSSESeqLenType::SHORT },
          { tk::kw::seq_med().string(), RNGSSESeqLenType::MEDIUM },
          { tk::kw::seq_long().string(), RNGSSESeqLenType::LONG } } ) {}
};

} // ctr::
} // tk::

#endif // RNGSSESeqLenOptions_h
