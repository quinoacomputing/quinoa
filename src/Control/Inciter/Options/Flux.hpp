// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Flux.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Flux function options for inciter
  \details   Flux function options for inciter
*/
// *****************************************************************************
#ifndef FluxOptions_h
#define FluxOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Flux types
enum class FluxType : uint8_t { LaxFriedrichs
                              , HLLC
                              , UPWIND
                              , AUSM
                              , HLL
                              , Rusanov
                              };

//! Pack/Unpack FluxType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, FluxType& e ) { PUP::pup( p, e ); }

//! \brief Flux options: outsource to base templated on enum type
class Flux : public tk::Toggle< FluxType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::laxfriedrichs
                                  , kw::hllc
                                  , kw::upwind
                                  , kw::ausm
                                  , kw::hll
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Flux() :
      tk::Toggle< FluxType >(
        //! Group, i.e., options, name
        kw::flux::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { FluxType::LaxFriedrichs, kw::laxfriedrichs::name() }
        , { FluxType::HLLC, kw::hllc::name() }
        , { FluxType::UPWIND, kw::upwind::name() }
        , { FluxType::AUSM, kw::ausm::name() }
        , { FluxType::HLL, kw::hll::name() }
        },
        //! keywords -> Enums
        { { kw::laxfriedrichs::string(), FluxType::LaxFriedrichs }
        , { kw::hllc::string(), FluxType::HLLC }
        , { kw::upwind::string(), FluxType::UPWIND }
        , { kw::ausm::string(), FluxType::AUSM }
        , { kw::hll::string(), FluxType::HLL }
        } )
    {}

};

} // ctr::
} // inciter::

#endif // FluxOptions_h
