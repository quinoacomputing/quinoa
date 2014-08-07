//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/CoeffPolicy.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:33:32 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
      Toggle< CoeffPolicyType >( "SDE coefficients policy",
        //! Enums -> names
        { { CoeffPolicyType::CONSTANT, kw::constant().name() } },
        //! keywords -> Enums
        {  { kw::constant().string(), CoeffPolicyType::CONSTANT } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaCoeffPolicyOptions_h
