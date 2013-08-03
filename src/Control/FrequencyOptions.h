//******************************************************************************
/*!
  \file      src/Control/FrequencyOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:41:44 2013
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

namespace Quinoa {

namespace select {

//! Frequency model types
enum class FrequencyType : uint8_t { NO_FREQUENCY=0,
                                     GAMMA };

//! Class with base templated on the above enum class with associations
class Frequency : public Toggle<FrequencyType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Frequency() : Toggle<FrequencyType>(names, values) {
      //! Enums -> names
      names[FrequencyType::NO_FREQUENCY] = "No frequency";
      names[FrequencyType::GAMMA] = "Gamma";
      //! keywords -> Enums
      values["no_frequency"] = FrequencyType::NO_FREQUENCY;
      values["freq_gamma"] = FrequencyType::GAMMA;
    }

  private:
    //! Don't permit copy constructor
    Frequency(const Frequency&) = delete;
    //! Don't permit copy assigment
    Frequency& operator=(const Frequency&) = delete;
    //! Don't permit move constructor
    Frequency(Frequency&&) = delete;
    //! Don't permit move assigment
    Frequency& operator=(Frequency&&) = delete;

    std::map<FrequencyType, std::string> names;
    std::map<std::string, FrequencyType> values;
};

} // namespace select

} // namespace Quinoa

#endif // FrequencyOptions_h
