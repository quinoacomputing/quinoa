//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/DiffEq.h
  \author    J. Bakosi
  \date      Thu 20 Nov 2014 09:22:23 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Differential equation options and associations
  \details   Differential equation options and associations
*/
//******************************************************************************
#ifndef QuinoaDiffEqOptions_h
#define QuinoaDiffEqOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>
#include <Quinoa/Options/InitPolicy.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {
namespace ctr {

//! Differential equation types
enum class DiffEqType : uint8_t { NO_DIFFEQ=0,
                                  OU,
                                  DIAG_OU,
                                  LOGNORMAL,
                                  SKEWNORMAL,
                                  GAMMA,
                                  BETA,
                                  DIRICHLET,
                                  GENDIR,
                                  WRIGHTFISHER };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, DiffEqType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using DiffEqKey =
  tk::tuple::tagged_tuple< tag::diffeq,      ctr::DiffEqType,
                           tag::initpolicy,  ctr::InitPolicyType,
                           tag::coeffpolicy, ctr::CoeffPolicyType >;

//! Class with base templated on the above enum class with associations
class DiffEq : public tk::Toggle< DiffEqType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit DiffEq() :
      Toggle< DiffEqType >( "Differential equation",
        //! Enums -> names
        { { DiffEqType::NO_DIFFEQ, "n/a" },
          { DiffEqType::OU, kw::ornstein_uhlenbeck().name() },
          { DiffEqType::DIAG_OU, kw::diag_ornstein_uhlenbeck().name() },
          { DiffEqType::LOGNORMAL, kw::lognormal().name() },
          { DiffEqType::SKEWNORMAL, kw::skewnormal().name() },
          { DiffEqType::GAMMA, kw::gamma().name() },
          { DiffEqType::BETA, kw::beta().name() },
          { DiffEqType::DIRICHLET, kw::dirichlet().name() },
          { DiffEqType::GENDIR, kw::gendir().name() },
          { DiffEqType::WRIGHTFISHER, kw::wrightfisher().name() } },
        //! keywords -> Enums
        { { "no_diffeq", DiffEqType::NO_DIFFEQ },
          { kw::ornstein_uhlenbeck().string(), DiffEqType::OU },
          { kw::diag_ornstein_uhlenbeck().string(), DiffEqType::DIAG_OU },
          { kw::lognormal().string(), DiffEqType::LOGNORMAL },
          { kw::skewnormal().string(), DiffEqType::SKEWNORMAL },
          { kw::gamma().string(), DiffEqType::GAMMA },
          { kw::beta().string(), DiffEqType::BETA },
          { kw::dirichlet().string(), DiffEqType::DIRICHLET },
          { kw::gendir().string(), DiffEqType::GENDIR },
          { kw::wrightfisher().string(), DiffEqType::WRIGHTFISHER } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaDiffEqOptions_h
