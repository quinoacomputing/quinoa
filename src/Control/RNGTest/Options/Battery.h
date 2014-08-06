//******************************************************************************
/*!
  \file      src/Control/RNGTest/Options/Battery.h
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 03:52:02 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator test batteries options and associations
  \details   Random number generator test batteries options and associations
*/
//******************************************************************************
#ifndef RNGTestBatteryOptions_h
#define RNGTestBatteryOptions_h

#include <map>

#include <Toggle.h>
#include <RNGTest/InputDeck/Keywords.h>
#include <PUPUtil.h>

namespace rngtest {
namespace ctr {

//! Random number generator battery types
enum class BatteryType : uint8_t { NO_BATTERY=0,
                                   SMALLCRUSH,
                                   CRUSH,
                                   BIGCRUSH };

//! Pack/Unpack BatteryType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, BatteryType& b ) { tk::pup( p, b ); }

//! Class with base templated on the above enum class with associations
class Battery : public tk::Toggle< BatteryType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Battery() :
      Toggle< BatteryType >( "RNG battery",
        //! Enums -> names
        { { BatteryType::NO_BATTERY, "n/a" },
          { BatteryType::SMALLCRUSH, kw::smallcrush().name() },
          { BatteryType::CRUSH, kw::crush().name() },
          { BatteryType::BIGCRUSH, kw::bigcrush().name() } },
        //! keywords -> Enums
        { { "no_battery", BatteryType::NO_BATTERY },
          { kw::smallcrush().string(), BatteryType::SMALLCRUSH },
          { kw::crush().string(), BatteryType::CRUSH },
          { kw::bigcrush().string(), BatteryType::BIGCRUSH } } ) {}
};

} // ctr::
} // rngtest::

#endif // RNGTestBatteryOptions_h
