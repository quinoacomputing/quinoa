//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/SDE.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 05:35:59 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     SDE options and associations
  \details   SDE options and associations
*/
//******************************************************************************
#ifndef QuinoaSDEOptions_h
#define QuinoaSDEOptions_h

#include <map>

#include <Model.h>
#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>
#include <Quinoa/Options/InitPolicy.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {
namespace ctr {

//! SDE types
enum class SDEType : uint8_t { NO_SDE=0,
                               ORNSTEIN_UHLENBECK,
                               LOGNORMAL,
                               SKEWNORMAL,
                               DIRICHLET,
                               GENDIR };

//! SDE factory type
using SDEKey = tk::tuple::tagged_tuple< tag::sde,         SDEType,
                                        tag::initpolicy,  InitPolicyType,
                                        tag::coeffpolicy, CoeffPolicyType >;
using SDEFactory = std::map< SDEKey, std::function< Model*() > >;

//! Class with base templated on the above enum class with associations
class SDE : public tk::Toggle< SDEType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit SDE() : Toggle< SDEType >( "SDE", names, values ) {}

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
    const kw::lognormal lognormal {};
    const kw::skewnormal skewnormal {};
    const kw::dirichlet dir {};
    const kw::gendir gendir {};

    //! Enums -> names
    const std::map<SDEType, std::string> names {
      { SDEType::NO_SDE, "n/a" },
      { SDEType::ORNSTEIN_UHLENBECK, ou.name() },
      { SDEType::LOGNORMAL, lognormal.name() },
      { SDEType::SKEWNORMAL, skewnormal.name() },
      { SDEType::DIRICHLET, dir.name() },
      { SDEType::GENDIR, gendir.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, SDEType> values {
      { "no_sde", SDEType::NO_SDE },
      { ou.string(), SDEType::ORNSTEIN_UHLENBECK },
      { skewnormal.string(), SDEType::SKEWNORMAL },
      { dir.string(), SDEType::DIRICHLET },
      { gendir.string(), SDEType::GENDIR }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaSDEOptions_h
