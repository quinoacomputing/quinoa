//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.h
  \author    J. Bakosi
  \date      Tue 14 Jan 2014 07:16:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mass model options and associations
  \details   Mass model options and associations
*/
//******************************************************************************
#ifndef QuinoaMassOptions_h
#define QuinoaMassOptions_h

#include <map>
#include <list>

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
class Mass : public tk::Toggle<MassType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Mass() : Toggle<MassType>("Mass", names, values) {}

  private:
    //! Don't permit copy constructor
    Mass(const Mass&) = delete;
    //! Don't permit copy assigment
    Mass& operator=(const Mass&) = delete;
    //! Don't permit move constructor
    Mass(Mass&&) = delete;
    //! Don't permit move assigment
    Mass& operator=(Mass&&) = delete;

    //! Get access to mass keywords
    const kw::mass_beta beta {};

    //! Enums -> names
    const std::map<MassType, std::string> names {
      { MassType::NO_MASS, "n/a" },
      { MassType::BETA, beta.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, MassType> values {
      { "no_mass", MassType::NO_MASS },
      { beta.string(), MassType::BETA }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaMassOptions_h
