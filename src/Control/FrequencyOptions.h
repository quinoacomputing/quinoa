//******************************************************************************
/*!
  \file      src/Control/FrequencyOptions.h
  \author    J. Bakosi
  \date      Wed May 29 07:22:08 2013
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
enum class FrequencyTypes { NO_FREQUENCY=0,
                            GAMMA };

//! Class with base templated on the above enum class with associations
class Frequency : public Toggle<FrequencyTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Frequency() : Toggle<FrequencyTypes>(names, values) {
      //! Enums -> names
      names[FrequencyTypes::NO_FREQUENCY] = "No energy";
      names[FrequencyTypes::GAMMA] = "Gamma";
      //! keywords -> Enums
      values["no_energy"] = FrequencyTypes::NO_FREQUENCY;
      values["gamma"] = FrequencyTypes::GAMMA;
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

    std::map<FrequencyTypes, std::string> names;
    std::map<std::string, FrequencyTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // FrequencyOptions_h
