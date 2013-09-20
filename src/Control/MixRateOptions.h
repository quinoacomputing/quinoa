//******************************************************************************
/*!
  \file      src/Control/MixRateOptions.h
  \author    J. Bakosi
  \date      Fri Sep 20 13:40:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef MixRateOptions_h
#define MixRateOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>
#include <QuinoaKeywords.h>

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
    explicit MixRate() : Toggle<MixRateType>(names, values) {}

  private:
    //! Don't permit copy constructor
    MixRate(const MixRate&) = delete;
    //! Don't permit copy assigment
    MixRate& operator=(const MixRate&) = delete;
    //! Don't permit move constructor
    MixRate(MixRate&&) = delete;
    //! Don't permit move assigment
    MixRate& operator=(MixRate&&) = delete;

    //! Enums -> names
    const std::map<MixRateType, std::string> names {
      { MixRateType::NO_MIXRATE, "" },
      { MixRateType::GAMMA, "Gamma" }
    };

    //! Get access to mixrate keywords
    const grm::kw::mixrate_gamma gamma {};

    //! keywords -> Enums
    const std::map<std::string, MixRateType> values {
      { "no_mixrate", MixRateType::NO_MIXRATE },
      { gamma.string(), MixRateType::GAMMA }
    };
};

} // sel::
} // quinoa::

#endif // MixRateOptions_h
