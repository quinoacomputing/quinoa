//******************************************************************************
/*!
  \file      src/Control/HydroOptions.h
  \author    J. Bakosi
  \date      Wed Sep 11 16:48:39 2013
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

namespace quinoa {

namespace select {

//! Hydro model types
enum class HydroType : uint8_t { NO_HYDRO=0,
                                 SLM,
                                 GLM };

//! Class with base templated on the above enum class with associations
class Hydro : public Toggle<HydroType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Hydro() : Toggle<HydroType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;

    //! Enums -> names
    const std::map<HydroType, std::string> names {
      { HydroType::NO_HYDRO, "" },
      { HydroType::SLM, "Simplified Langevin" },
      { HydroType::GLM, "Generalized Langevin"}
    };

    //! keywords -> Enums
    const std::map<std::string, HydroType> values {
      { "no_hydro", HydroType::NO_HYDRO },
      { "hydro_slm", HydroType::SLM },
      { "hydro_glm", HydroType::GLM }
    };
};

} // namespace select

} // namespace quinoa

#endif // HydroOptions_h
