//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:41:04 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Mass model options and associations
  \details   Mass model options and associations
*/
//******************************************************************************
#ifndef QuinoaMassOptions_h
#define QuinoaMassOptions_h

#include <map>

#include <Model.h>
#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Mass model types
enum class MassType : uint8_t { NO_MASS=0,
                                BETA };

//! Mass model factory type
using MassFactory = std::map< MassType, std::function<Model*()> >;

//! Class with base templated on the above enum class with associations
class Mass : public tk::Toggle< MassType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Mass() :
      Toggle< MassType >( "Mass",
        //! Enums -> names
        { { MassType::NO_MASS, "n/a" },
          { MassType::BETA, kw::mass_beta().name() } },
        //! keywords -> Enums
        { { "no_mass", MassType::NO_MASS },
          { kw::mass_beta().string(), MassType::BETA } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaMassOptions_h
