// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Limiter.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Limiter options for DG
  \details   Limiter options for DG
*/
// *****************************************************************************
#ifndef LimiterOptions_h
#define LimiterOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Limiter types
enum class LimiterType : uint8_t { NOLIMITER
                                 , WENOP1 };

//! Pack/Unpack LimiterType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, LimiterType& e ) { PUP::pup( p, e ); }

//! \brief Limiter options: outsource to base templated on enum type
class Limiter : public tk::Toggle< LimiterType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::nolimiter
                                  , kw::wenop1
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Limiter() :
      tk::Toggle< LimiterType >(
        //! Group, i.e., options, name
        kw::limiter::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { LimiterType::NOLIMITER, kw::nolimiter::name() },
          { LimiterType::WENOP1, kw::wenop1::name() } },
        //! keywords -> Enums
        { { kw::nolimiter::string(), LimiterType::NOLIMITER },
          { kw::wenop1::string(), LimiterType::WENOP1 } } )
    {}

};

} // ctr::
} // inciter::

#endif // LimiterOptions_h
