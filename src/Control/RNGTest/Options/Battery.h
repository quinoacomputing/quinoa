//******************************************************************************
/*!
  \file      src/Control/RNGTest/Options/Battery.h
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 12:39:54 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test batteries options and associations
  \details   Random number generator test batteries options and associations
*/
//******************************************************************************
#ifndef RNGTestBatteryOptions_h
#define RNGTestBatteryOptions_h

#include <map>
#include <list>

#include <Toggle.h>
#include <Battery.h>

namespace rngtest {
namespace ctr {

//! Random number generator battery types
enum class BatteryType : uint8_t { NO_BATTERY=0,
                                   SMALLCRUSH,
                                   CRUSH,
                                   BIGCRUSH };

//! Battery factory type
using BatteryFactory = std::map< BatteryType, std::function<Battery*()> >;

//! Class with base templated on the above enum class with associations
class Battery : public tk::Toggle<BatteryType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Battery() :
      Toggle<BatteryType>
        ("RNG battery", names, values) {}

    //! Register batteries into factory
    void initFactory( BatteryFactory& f, std::list<std::string>& reg ) const;

  private:
    //! Don't permit copy constructor
    Battery(const Battery&) = delete;
    //! Don't permit copy assigment
    Battery& operator=(const Battery&) = delete;
    //! Don't permit move constructor
    Battery(Battery&&) = delete;
    //! Don't permit move assigment
    Battery& operator=(Battery&&) = delete;

    //! Enums -> names
    const std::map<BatteryType, std::string> names {
      { BatteryType::NO_BATTERY, "No battery" },
      { BatteryType::SMALLCRUSH, "SmallCrush" },
      { BatteryType::CRUSH, "Crush" },
      { BatteryType::BIGCRUSH, "BigCrush" }
    };

    //! keywords -> Enums
    const std::map<std::string, BatteryType> values {
      { "no_battery", BatteryType::NO_BATTERY },
      { "smallcrush", BatteryType::SMALLCRUSH },
      { "crush", BatteryType::CRUSH },
      { "bigcrush", BatteryType::BIGCRUSH }
    };
};

} // ctr::
} // rngtest::

#endif // RNGTestBatteryOptions_h
