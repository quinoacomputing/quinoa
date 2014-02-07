//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/CoeffPolicy.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 05:35:20 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE coefficients policy options and associations
  \details   SDE coefficients policy options and associations
*/
//******************************************************************************
#ifndef QuinoaCoeffPolicyOptions_h
#define QuinoaCoeffPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! SDE coefficients policy types
enum class CoeffPolicyType : uint8_t { CONSTANT=0 };

//! Class with base templated on the above enum class with associations
class CoeffPolicy : public tk::Toggle< CoeffPolicyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit CoeffPolicy() :
      Toggle< CoeffPolicyType >
            ( "SDE coefficients policy", names, values ) {}

  private:
    //! Don't permit copy constructor
    CoeffPolicy(const CoeffPolicy&) = delete;
    //! Don't permit copy assigment
    CoeffPolicy& operator=(const CoeffPolicy&) = delete;
    //! Don't permit move constructor
    CoeffPolicy(CoeffPolicy&&) = delete;
    //! Don't permit move assigment
    CoeffPolicy& operator=(CoeffPolicy&&) = delete;

    //! Get access to coefficients policy keywords
    const kw::constant constant {};

    //! Enums -> names
    const std::map<CoeffPolicyType, std::string> names {
      { CoeffPolicyType::CONSTANT, constant.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, CoeffPolicyType> values {
      { constant.string(), CoeffPolicyType::CONSTANT }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaCoeffPolicyOptions_h
