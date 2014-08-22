//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Mon 11 Aug 2014 09:34:39 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Differential equation initialization policy options
  \details   Differential equation initialization policy options
*/
//******************************************************************************
#ifndef QuinoaInitPolicyOptions_h
#define QuinoaInitPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Differential equation initializion policies
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
