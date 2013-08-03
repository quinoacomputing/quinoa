//******************************************************************************
/*!
  \file      src/Control/HydroOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:40:37 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model options and associations
  \details   Hydro model options and associations
*/
//******************************************************************************
#ifndef HydroOptions_h
#define HydroOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Hydro model types
enum class HydroType : uint8_t { NO_HYDRO=0,
                                 SLM,
                                 GLM };

//! Class with base templated on the above enum class with associations
class Hydro : public Toggle<HydroType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Hydro() : Toggle<HydroType>(names, values) {
      //! Enums -> names
      names[HydroType::NO_HYDRO] = "No hydro";
      names[HydroType::SLM] = "Simplified Langevin";
      names[HydroType::GLM] = "Generalized Langevin";
      //! keywords -> Enums
      values["no_hydro"] = HydroType::NO_HYDRO;
      values["hydro_slm"] = HydroType::SLM;
      values["hydro_glm"] = HydroType::GLM;
    }

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;

    std::map<HydroType, std::string> names;
    std::map<std::string, HydroType> values;
};

} // namespace select

} // namespace Quinoa

#endif // HydroOptions_h
