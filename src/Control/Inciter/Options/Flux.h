// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Flux.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Flux function options for inciter
  \details   Flux function options for inciter
*/
// *****************************************************************************
#ifndef FluxOptions_h
#define FluxOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Flux types
enum class FluxType : uint8_t { LaxFriedrichs
                              , HLLC };

//! Pack/Unpack FluxType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, FluxType& e ) { PUP::pup( p, e ); }

//! \brief Flux options: outsource to base templated on enum type
class Flux : public tk::Toggle< FluxType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::laxfriedrichs
                                  , kw::hllc
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Flux() :
      tk::Toggle< FluxType >(
        //! Group, i.e., options, name
        kw::flux::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { FluxType::LaxFriedrichs, kw::laxfriedrichs::name() },
          { FluxType::HLLC, kw::hllc::name() } },
        //! keywords -> Enums
        { { kw::laxfriedrichs::string(), FluxType::LaxFriedrichs },
          { kw::hllc::string(), FluxType::HLLC } } )
    {}

};

} // ctr::
} // inciter::

#endif // FluxOptions_h
