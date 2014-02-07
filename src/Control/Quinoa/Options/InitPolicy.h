//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 04:33:49 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
      Toggle< InitPolicyType >( "Initialization Policy", names, values ) {}

  private:
    //! Don't permit copy constructor
    InitPolicy(const InitPolicy&) = delete;
    //! Don't permit copy assigment
    InitPolicy& operator=(const InitPolicy&) = delete;
    //! Don't permit move constructor
    InitPolicy(InitPolicy&&) = delete;
    //! Don't permit move assigment
    InitPolicy& operator=(InitPolicy&&) = delete;

    //! Get access to InitPolicy keywords
    const kw::raw raw {};
    const kw::zero zero {};

    //! Enums -> names
    const std::map<InitPolicyType, std::string> names {
      { InitPolicyType::RAW, raw.name() },
      { InitPolicyType::ZERO, zero.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, InitPolicyType> values {
      { raw.string(), InitPolicyType::RAW },
      { zero.string(), InitPolicyType::ZERO }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaInitPolicyOptions_h
