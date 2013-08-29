//******************************************************************************
/*!
  \file      src/Control/FrequencyOptions.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:03:58 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef FrequencyOptions_h
#define FrequencyOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace quinoa {

namespace select {

//! Frequency model types
enum class FrequencyType : uint8_t { NO_FREQUENCY=0,
                                     GAMMA };

//! Class with base templated on the above enum class with associations
class Frequency : public Toggle<FrequencyType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Frequency() : Toggle<FrequencyType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Frequency(const Frequency&) = delete;
    //! Don't permit copy assigment
    Frequency& operator=(const Frequency&) = delete;
    //! Don't permit move constructor
    Frequency(Frequency&&) = delete;
    //! Don't permit move assigment
    Frequency& operator=(Frequency&&) = delete;

    //! Enums -> names
    const std::map<FrequencyType, std::string> names {
      { FrequencyType::NO_FREQUENCY, "No frequency" },
      { FrequencyType::GAMMA, "Gamma" }
    };

    //! keywords -> Enums
    const std::map<std::string, FrequencyType> values {
      { "no_frequency", FrequencyType::NO_FREQUENCY },
      { "freq_gamma", FrequencyType::GAMMA }
    };
};

} // namespace select

} // namespace quinoa

#endif // FrequencyOptions_h
