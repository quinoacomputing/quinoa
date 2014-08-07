//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/SDE.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:28:40 PM MDT
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
    explicit SDE() :
      Toggle< SDEType >( "SDE",
        //! Enums -> names
        { { SDEType::NO_SDE, "n/a" },
          { SDEType::ORNSTEIN_UHLENBECK, kw::ornstein_uhlenbeck().name() },
          { SDEType::LOGNORMAL, kw::lognormal().name() },
          { SDEType::SKEWNORMAL, kw::skewnormal().name() },
          { SDEType::DIRICHLET, kw::dirichlet().name() },
          { SDEType::GENDIR, kw::gendir().name() } },
        //! keywords -> Enums
        { { "no_sde", SDEType::NO_SDE },
          { kw::ornstein_uhlenbeck().string(), SDEType::ORNSTEIN_UHLENBECK },
          { kw::lognormal().string(), SDEType::LOGNORMAL },
          { kw::skewnormal().string(), SDEType::SKEWNORMAL },
          { kw::dirichlet().string(), SDEType::DIRICHLET },
          { kw::gendir().string(), SDEType::GENDIR } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaSDEOptions_h
