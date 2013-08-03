//******************************************************************************
/*!
  \file      src/Control/MixRateOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:41:10 2013
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
enum class MixRateType : uint8_t { NO_MIXRATE=0,
                                   GAMMA };

//! Class with base templated on the above enum class with associations
class MixRate : public Toggle<MixRateType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    MixRate() : Toggle<MixRateType>(names, values) {
      //! Enums -> names
      names[MixRateType::NO_MIXRATE] = "No mix rate";
      names[MixRateType::GAMMA] = "Gamma";
      //! keywords -> Enums
      values["no_mixrate"] = MixRateType::NO_MIXRATE;
      values["gamma"] = MixRateType::GAMMA;
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

    std::map<MixRateType, std::string> names;
    std::map<std::string, MixRateType> values;
};

} // namespace select

} // namespace Quinoa

#endif // MixRateOptions_h
