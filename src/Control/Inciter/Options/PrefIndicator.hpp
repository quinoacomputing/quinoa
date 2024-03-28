// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/PrefIndicator.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Options for adaptive indicators for p-adaptive DG scheme.
  \details   Options for adaptive indicators for p-adaptive DG scheme.
*/
// *****************************************************************************
#ifndef InciterPrefIndicatorOptions_h
#define InciterPrefIndicatorOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Types of adaptive indicators
enum class PrefIndicatorType : uint8_t { SPECTRAL_DECAY
                                       , NON_CONFORMITY };

//! Pack/Unpack PrefIndicatorType: forward overload to generic enum class
//! packer
inline void operator|( PUP::er& p, PrefIndicatorType& e ) { PUP::pup( p, e ); }

//! PrefIndicator options: outsource searches to base templated on enum type
class PrefIndicator : public tk::Toggle< PrefIndicatorType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PrefIndicator() :
      tk::Toggle< PrefIndicatorType >(
        //! Group, i.e., options, name
        "p-adaptive indicator",
        //! Enums -> names
        { { PrefIndicatorType::SPECTRAL_DECAY, "spectral_decay" },
          { PrefIndicatorType::NON_CONFORMITY, "non_conformity" } },
        //! keywords -> Enums
        { { "spectral_decay", PrefIndicatorType::SPECTRAL_DECAY },
          { "non_conformity", PrefIndicatorType::NON_CONFORMITY } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterPrefIndicatorOptions_h
