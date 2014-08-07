//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:39:45 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     SDE initialization policy options and associations
  \details   SDE initialization policy options and associations
*/
//******************************************************************************
#ifndef QuinoaInitPolicyOptions_h
#define QuinoaInitPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! InitPolicy types
enum class InitPolicyType : uint8_t { RAW=0,
                                      ZERO };

//! Class with base templated on the above enum class with associations
class InitPolicy : public tk::Toggle< InitPolicyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit InitPolicy() :
      Toggle< InitPolicyType >( "Initialization Policy",
        //! Enums -> names
        { { InitPolicyType::RAW, kw::raw().name() },
          { InitPolicyType::ZERO, kw::zero().name() } },
        //! keywords -> Enums
        { { kw::raw().string(), InitPolicyType::RAW },
          { kw::zero().string(), InitPolicyType::ZERO } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaInitPolicyOptions_h
