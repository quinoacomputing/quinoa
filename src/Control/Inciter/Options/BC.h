// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/BC.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Boundary condition options for inciter
  \details   Boundary condition options for inciter
*/
// *****************************************************************************
#ifndef InciterBCOptions_h
#define InciterBCOptions_h

#include "NoWarning/vector.h"

#include "TaggedTuple.h"
#include "Toggle.h"
#include "Keywords.h"

namespace inciter {
namespace ctr {

//! Boundary condition types
enum class BCType : uint8_t { SYM=0,
                              INLET,
                              OUTLET };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, BCType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class BC : public tk::Toggle< BCType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::bc_sym
                                       , kw::bc_inlet
                                       , kw::bc_outlet
                                       >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit BC() :
      tk::Toggle< BCType >( "Partial differential equation",
        //! Enums -> names
        { { BCType::SYM, kw::bc_sym::name() },
          { BCType::INLET, kw::bc_inlet::name() },
          { BCType::OUTLET, kw::bc_outlet::name() } },
        //! keywords -> Enums
        { { kw::bc_sym::string(), BCType::SYM },
          { kw::bc_inlet::string(), BCType::INLET },
          { kw::bc_outlet::string(), BCType::OUTLET } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterBCOptions_h
