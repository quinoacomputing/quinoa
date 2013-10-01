//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MixRate.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:08:46 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef QuinoaMixRateOptions_h
#define QuinoaMixRateOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace sel {

//! Material mix rate model types
enum class MixRateType : uint8_t { NO_MIXRATE=0,
                                   GAMMA };

//! Class with base templated on the above enum class with associations
class MixRate : public Toggle<MixRateType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MixRate() :
      Toggle<MixRateType>("Material mix rate", names, values) {}

  private:
    //! Don't permit copy constructor
    MixRate(const MixRate&) = delete;
    //! Don't permit copy assigment
    MixRate& operator=(const MixRate&) = delete;
    //! Don't permit move constructor
    MixRate(MixRate&&) = delete;
    //! Don't permit move assigment
    MixRate& operator=(MixRate&&) = delete;

    //! Get access to mixrate keywords
    const grm::kw::mixrate_gamma gamma {};

    //! Enums -> names
    const std::map<MixRateType, std::string> names {
      { MixRateType::NO_MIXRATE, "n/a" },
      { MixRateType::GAMMA, gamma.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, MixRateType> values {
      { "no_mixrate", MixRateType::NO_MIXRATE },
      { gamma.string(), MixRateType::GAMMA }
    };
};

} // sel::
} // quinoa::

#endif // QuinoaMixRateOptions_h
