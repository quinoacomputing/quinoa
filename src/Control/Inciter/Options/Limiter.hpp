// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Limiter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Limiter options for DG
  \details   Limiter options for DG
*/
// *****************************************************************************
#ifndef LimiterOptions_h
#define LimiterOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Limiter types
enum class LimiterType : uint8_t { NOLIMITER
                                 , WENOP1
                                 , SUPERBEEP1
                                 , VERTEXBASEDP1 };

//! Pack/Unpack LimiterType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, LimiterType& e ) { PUP::pup( p, e ); }

//! \brief Limiter options: outsource to base templated on enum type
class Limiter : public tk::Toggle< LimiterType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Limiter() :
      tk::Toggle< LimiterType >(
        //! Group, i.e., options, name
        "Limiter",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { LimiterType::NOLIMITER, "nolimiter" },
          { LimiterType::WENOP1, "wenop1" },
          { LimiterType::SUPERBEEP1, "superbeep1" },
          { LimiterType::VERTEXBASEDP1, "vertexbasedp1" } },
        //! keywords -> Enums
        { { "nolimiter", LimiterType::NOLIMITER },
          { "wenop1", LimiterType::WENOP1 },
          { "superbeep1", LimiterType::SUPERBEEP1 },
          { "vertexbasedp1", LimiterType::VERTEXBASEDP1 } } )
    {}

};

} // ctr::
} // inciter::

#endif // LimiterOptions_h
