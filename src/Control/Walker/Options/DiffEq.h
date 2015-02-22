//******************************************************************************
/*!
  \file      src/Control/Walker/Options/DiffEq.h
  \author    J. Bakosi
  \date      Mon 09 Feb 2015 12:33:53 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Differential equation options and associations
  \details   Differential equation options and associations
*/
//******************************************************************************
#ifndef WalkerDiffEqOptions_h
#define WalkerDiffEqOptions_h

#include <boost/mpl/vector.hpp>

#include <TaggedTuple.h>
#include <Toggle.h>
#include <Keywords.h>
#include <Options/InitPolicy.h>
#include <Options/CoeffPolicy.h>

namespace walker {
namespace ctr {

//! Differential equation types
enum class DiffEqType : uint8_t { NO_DIFFEQ=0,
                                  OU,
                                  DIAG_OU,
                                  SKEWNORMAL,
                                  GAMMA,
                                  BETA,
                                  FUNCBETA,
                                  DIRICHLET,
                                  GENDIR,
                                  WRIGHTFISHER };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, DiffEqType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using DiffEqKey =
  tk::tuple::tagged_tuple< tag::diffeq,      ctr::DiffEqType,
                           tag::initpolicy,  tk::ctr::InitPolicyType,
                           tag::coeffpolicy, tk::ctr::CoeffPolicyType >;

//! Class with base templated on the above enum class with associations
class DiffEq : public tk::Toggle< DiffEqType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::ornstein_uhlenbeck
                                       , kw::diag_ou
                                       , kw::skewnormal
                                       , kw::gamma
                                       , kw::beta
                                       , kw::funcbeta
                                       , kw::dirichlet
                                       , kw::gendir
                                       , kw::wrightfisher
                                       >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit DiffEq() :
      Toggle< DiffEqType >( "Differential equation",
        //! Enums -> names
        { { DiffEqType::NO_DIFFEQ, "n/a" },
          { DiffEqType::OU, kw::ornstein_uhlenbeck::name() },
          { DiffEqType::DIAG_OU, kw::diag_ou::name() },
          { DiffEqType::SKEWNORMAL, kw::skewnormal::name() },
          { DiffEqType::GAMMA, kw::gamma::name() },
          { DiffEqType::BETA, kw::beta::name() },
          { DiffEqType::FUNCBETA, kw::funcbeta::name() },
          { DiffEqType::DIRICHLET, kw::dirichlet::name() },
          { DiffEqType::GENDIR, kw::gendir::name() },
          { DiffEqType::WRIGHTFISHER, kw::wrightfisher::name() } },
        //! keywords -> Enums
        { { "no_diffeq", DiffEqType::NO_DIFFEQ },
          { kw::ornstein_uhlenbeck::string(), DiffEqType::OU },
          { kw::diag_ou::string(), DiffEqType::DIAG_OU },
          { kw::skewnormal::string(), DiffEqType::SKEWNORMAL },
          { kw::gamma::string(), DiffEqType::GAMMA },
          { kw::beta::string(), DiffEqType::BETA },
          { kw::funcbeta::string(), DiffEqType::FUNCBETA },
          { kw::dirichlet::string(), DiffEqType::DIRICHLET },
          { kw::gendir::string(), DiffEqType::GENDIR },
          { kw::wrightfisher::string(), DiffEqType::WRIGHTFISHER } } ) {}
};

} // ctr::
} // walker::

#endif // WalkerDiffEqOptions_h
