//******************************************************************************
/*!
  \file      src/Control/Inciter/Options/DiffEq.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:14:59 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Differential equation options and associations for inciter
  \details   Differential equation options and associations for inciter
*/
//******************************************************************************
#ifndef InciterDiffEqOptions_h
#define InciterDiffEqOptions_h

#include <boost/mpl/vector.hpp>

#include "TaggedTuple.h"
#include "Toggle.h"
#include "Keywords.h"
#include "Inciter/Options/Problem.h"

namespace inciter {
namespace ctr {

//! Differential equation types
enum class DiffEqType : uint8_t { NO_DIFFEQ=0,
                                  SCALAR };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, DiffEqType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using DiffEqKey =
  tk::tuple::tagged_tuple< tag::diffeq,  DiffEqType,
                           tag::problem, ProblemType >;

//! Class with base templated on the above enum class with associations
class DiffEq : public tk::Toggle< DiffEqType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::scalar
                                       >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit DiffEq() :
      Toggle< DiffEqType >( "Differential equation",
        //! Enums -> names
        { { DiffEqType::NO_DIFFEQ, "n/a" },
          { DiffEqType::SCALAR, kw::scalar::name() } },
        //! keywords -> Enums
        { { "no_diffeq", DiffEqType::NO_DIFFEQ },
          { kw::scalar::string(), DiffEqType::SCALAR } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterDiffEqOptions_h
