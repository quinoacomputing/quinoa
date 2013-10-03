//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.h
  \author    J. Bakosi
  \date      Thu Oct  3 17:39:26 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model options and associations
  \details   Hydro model options and associations
*/
//******************************************************************************
#ifndef QuinoaHydroOptions_h
#define QuinoaHydroOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace sel {

//! Hydro model types
enum class HydroType : uint8_t { NO_HYDRO=0,
                                 SLM,
                                 GLM };

//! Class with base templated on the above enum class with associations
class Hydro : public Toggle<HydroType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Hydro() : Toggle<HydroType>("Hydrodynamics", names, values) {}

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;

    //! Get access to hydro keywords
    const kw::hydro_slm slm {};
    const kw::hydro_glm glm {};

    //! Enums -> names
    const std::map<HydroType, std::string> names {
      { HydroType::NO_HYDRO, "n/a" },
      { HydroType::SLM, slm.name() },
      { HydroType::GLM, glm.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, HydroType> values {
      { "no_hydro", HydroType::NO_HYDRO },
      { slm.string(), HydroType::SLM },
      { glm.string(), HydroType::GLM }
    };
};

} // sel::
} // quinoa::

#endif // QuinoaHydroOptions_h
