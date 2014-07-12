//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 04:32:35 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Hydro model options and associations
  \details   Hydro model options and associations
*/
//******************************************************************************
#ifndef QuinoaHydroOptions_h
#define QuinoaHydroOptions_h

#include <map>

#include <Model.h>
#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Hydro model types
enum class HydroType : uint8_t { NO_HYDRO=0,
                                 SLM,
                                 GLM };

//! Hydrodynamics model factory type
using HydroFactory = std::map< HydroType, std::function<Model*()> >;

//! Class with base templated on the above enum class with associations
class Hydro : public tk::Toggle<HydroType> {

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

} // ctr::
} // quinoa::

#endif // QuinoaHydroOptions_h
