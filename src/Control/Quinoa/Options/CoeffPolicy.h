//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/CoeffPolicy.h
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 10:09:09 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Differential equation coefficients policy options and associations
  \details   Differential equation coefficients policy options and associations
*/
//******************************************************************************
#ifndef QuinoaCoeffPolicyOptions_h
#define QuinoaCoeffPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Differential equation coefficients policies
enum class CoeffPolicyType : uint8_t { CONSTANT=0 };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, CoeffPolicyType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class CoeffPolicy : public tk::Toggle< CoeffPolicyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit CoeffPolicy() :
      Toggle< CoeffPolicyType >( "Coefficients Policy",
        //! Enums -> names
        { { CoeffPolicyType::CONSTANT, kw::constant().name() } },
        //! keywords -> Enums
        {  { kw::constant().string(), CoeffPolicyType::CONSTANT } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaCoeffPolicyOptions_h
