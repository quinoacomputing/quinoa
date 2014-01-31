//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/SDE.h
  \author    J. Bakosi
  \date      Fri Jan 31 09:49:01 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE options and associations
  \details   SDE options and associations
*/
//******************************************************************************
#ifndef QuinoaSDEOptions_h
#define QuinoaSDEOptions_h

#include <map>
#include <list>

#include <Model.h>
#include <Toggle.h>

namespace quinoa {
namespace ctr {

//! SDE types
enum class SDEType : uint8_t { NO_SDE=0,
                               ORNSTEIN_UHLENBECK,
                               DIRICHLET,
                               GENDIR };

//! SDE factory type
using SDEFactory = std::map< SDEType, std::function<Model*()> >;

//! Class with base templated on the above enum class with associations
class SDE : public tk::Toggle<SDEType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit SDE() : Toggle<SDEType>("SDE", names, values) {}

  private:
    //! Don't permit copy constructor
    SDE(const SDE&) = delete;
    //! Don't permit copy assigment
    SDE& operator=(const SDE&) = delete;
    //! Don't permit move constructor
    SDE(SDE&&) = delete;
    //! Don't permit move assigment
    SDE& operator=(SDE&&) = delete;

    //! Get access to SDE keywords
    const kw::ornstein_uhlenbeck ou {};
    const kw::dirichlet dir {};
    const kw::gendir gendir {};

    //! Enums -> names
    const std::map<SDEType, std::string> names {
      { SDEType::NO_SDE, "n/a" },
      { SDEType::ORNSTEIN_UHLENBECK, ou.name() },
      { SDEType::DIRICHLET, dir.name() },
      { SDEType::GENDIR, gendir.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, SDEType> values {
      { "no_sde", SDEType::NO_SDE },
      { ou.string(), SDEType::ORNSTEIN_UHLENBECK },
      { dir.string(), SDEType::DIRICHLET },
      { gendir.string(), SDEType::GENDIR }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaSDEOptions_h
