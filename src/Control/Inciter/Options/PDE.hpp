// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/PDE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Partial differential equation options and associations for inciter
  \details   Partial differential equation options and associations for inciter
*/
// *****************************************************************************
#ifndef InciterPDEOptions_h
#define InciterPDEOptions_h

#include <brigand/sequences/list.hpp>

#include "TaggedTuple.hpp"
#include "Toggle.hpp"
#include "Keywords.hpp"
#include "Inciter/Options/Physics.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {
namespace ctr {

//! Differential equation types
enum class PDEType : uint8_t { TRANSPORT,
                               COMPFLOW,
                               MULTIMAT };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PDEType& e ) { PUP::pup( p, e ); }

//! Differential equation key used to access a diff eq in a factory
using PDEKey =
  tk::TaggedTuple< brigand::list<
      tag::pde,         PDEType
    , tag::physics,     ctr::PhysicsType
    , tag::problem,     ctr::ProblemType
  > >;

//! Class with base templated on the above enum class with associations
class PDE : public tk::Toggle< PDEType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::transport
                                  , kw::compflow
                                  , kw::multimat
                                  >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit PDE() :
      tk::Toggle< PDEType >( "Partial differential equation",
        //! Enums -> names
        { { PDEType::TRANSPORT, kw::transport::name() },
          { PDEType::COMPFLOW, kw::compflow::name() },
          { PDEType::MULTIMAT, kw::multimat::name() } },
        //! keywords -> Enums
        { { kw::transport::string(), PDEType::TRANSPORT },
          { kw::compflow::string(), PDEType::COMPFLOW },
          { kw::multimat::string(), PDEType::MULTIMAT } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterPDEOptions_h
