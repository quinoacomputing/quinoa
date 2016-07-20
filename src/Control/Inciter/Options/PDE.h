// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/PDE.h
  \author    J. Bakosi
  \date      Tue Jul 19 22:22:55 2016
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Partial differential equation options and associations for inciter
  \details   Partial differential equation options and associations for inciter
*/
// *****************************************************************************
#ifndef InciterPDEOptions_h
#define InciterPDEOptions_h

#include <boost/mpl/vector.hpp>

#include "TaggedTuple.h"
#include "Toggle.h"
#include "Keywords.h"
#include "Inciter/Options/Problem.h"

namespace inciter {
namespace ctr {

//! Differential equation types
enum class PDEType : uint8_t { NO_PDE=0,
                               ADV_DIFF,
                               EULER,
                               COMPNS };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDEType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using PDEKey =
  tk::tuple::tagged_tuple< tag::pde,         PDEType,
                           tag::problem,     ctr::ProblemType >;

//! Class with base templated on the above enum class with associations
class PDE : public tk::Toggle< PDEType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::advdiff
                                       , kw::euler
                                       , kw::compns
                                       >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDE() :
      Toggle< PDEType >( "Partial differential equation",
        //! Enums -> names
        { { PDEType::NO_PDE, "n/a" },
          { PDEType::ADV_DIFF, kw::advdiff::name() },
          { PDEType::EULER, kw::euler::name() },
          { PDEType::COMPNS, kw::compns::name() } },
        //! keywords -> Enums
        { { "no_pde", PDEType::NO_PDE },
          { kw::advdiff::string(), PDEType::ADV_DIFF },
          { kw::euler::string(), PDEType::EULER },
          { kw::compns::string(), PDEType::COMPNS } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterPDEOptions_h
