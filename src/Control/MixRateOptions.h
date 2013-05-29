//******************************************************************************
/*!
  \file      src/Control/MixRateOptions.h
  \author    J. Bakosi
  \date      Wed May 29 07:27:48 2013
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

namespace Quinoa {

namespace select {

//! Material mix rate model types
enum class MixRateTypes { NO_MIXRATE=0,
                          GAMMA };

//! Class with base templated on the above enum class with associations
class MixRate : public Toggle<MixRateTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    MixRate() : Toggle<MixRateTypes>(names, values) {
      //! Enums -> names
      names[MixRateTypes::NO_MIXRATE] = "No mix rate";
      names[MixRateTypes::GAMMA] = "Gamma";
      //! keywords -> Enums
      values["no_mixrate"] = MixRateTypes::NO_MIXRATE;
      values["gamma"] = MixRateTypes::GAMMA;
    }

  private:
    //! Don't permit copy constructor
    MixRate(const MixRate&) = delete;
    //! Don't permit copy assigment
    MixRate& operator=(const MixRate&) = delete;
    //! Don't permit move constructor
    MixRate(MixRate&&) = delete;
    //! Don't permit move assigment
    MixRate& operator=(MixRate&&) = delete;

    std::map<MixRateTypes, std::string> names;
    std::map<std::string, MixRateTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // MixRateOptions_h
