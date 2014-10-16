//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.h
  \author    J. Bakosi
  \date      Tue 16 Sep 2014 08:14:16 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
namespace ctr {

//! Hydro model types
enum class HydroType : uint8_t { NO_HYDRO=0,
                                 SLM,
                                 GLM };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, HydroType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class Hydro : public tk::Toggle< HydroType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Hydro() :
      Toggle< HydroType >( "Hydrodynamics",
        //! Enums -> names
        { { HydroType::NO_HYDRO, "n/a" },
          { HydroType::SLM, kw::hydro_slm().name() },
          { HydroType::GLM, kw::hydro_glm().name() } },
        //! keywords -> Enums
        { { "no_hydro", HydroType::NO_HYDRO },
          { kw::hydro_slm().string(), HydroType::SLM },
          { kw::hydro_glm().string(), HydroType::GLM } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaHydroOptions_h
