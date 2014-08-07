//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MixRate.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:44:41 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef QuinoaMixRateOptions_h
#define QuinoaMixRateOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Material mix rate model types
enum class MixRateType : uint8_t { NO_MIXRATE=0,
                                   GAMMA };

//! Class with base templated on the above enum class with associations
class MixRate : public tk::Toggle<MixRateType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MixRate() :
      Toggle< MixRateType >( "Material mix rate",
      //! Enums -> names
      { { MixRateType::NO_MIXRATE, "n/a" },
        { MixRateType::GAMMA, kw::mixrate_gamma().name() } },
      //! keywords -> Enums
      { { "no_mixrate", MixRateType::NO_MIXRATE },
        { kw::mixrate_gamma().string(), MixRateType::GAMMA } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaMixRateOptions_h
