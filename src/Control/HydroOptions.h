//******************************************************************************
/*!
  \file      src/Control/HydroOptions.h
  \author    J. Bakosi
  \date      Fri May 31 13:29:50 2013
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
enum class HydroTypes : uint8_t { NO_HYDRO=0,
                                  SLM,
                                  GLM };

//! Class with base templated on the above enum class with associations
class Hydro : public Toggle<HydroTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Hydro() : Toggle<HydroTypes>(names, values) {
      //! Enums -> names
      names[HydroTypes::NO_HYDRO] = "No hydro";
      names[HydroTypes::SLM] = "Simplified Langevin";
      names[HydroTypes::GLM] = "Generalized Langevin";
      //! keywords -> Enums
      values["no_hydro"] = HydroTypes::NO_HYDRO;
      values["hydro_slm"] = HydroTypes::SLM;
      values["hydro_glm"] = HydroTypes::GLM;
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

    std::map<HydroTypes, std::string> names;
    std::map<std::string, HydroTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // HydroOptions_h
