// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/PDE.h
  \author    J. Bakosi
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
#include "Inciter/Options/Physics.h"
#include "Inciter/Options/Problem.h"

namespace inciter {
namespace ctr {

//! Differential equation types
enum class PDEType : uint8_t { TRANSPORT=0,
                               COMPFLOW };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDEType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using PDEKey =
  tk::tuple::tagged_tuple< tag::pde,         PDEType,
                           tag::physics,     ctr::PhysicsType,
                           tag::problem,     ctr::ProblemType >;

//! Class with base templated on the above enum class with associations
class PDE : public tk::Toggle< PDEType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::transport
                                       , kw::poisson
                                       , kw::compflow
                                       >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDE() :
      tk::Toggle< PDEType >( "Partial differential equation",
        //! Enums -> names
        { { PDEType::TRANSPORT, kw::transport::name() },
          { PDEType::COMPFLOW, kw::compflow::name() } },
        //! keywords -> Enums
        { { kw::transport::string(), PDEType::TRANSPORT },
          { kw::compflow::string(), PDEType::COMPFLOW } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterPDEOptions_h
