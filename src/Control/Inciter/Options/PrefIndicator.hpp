// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/PrefIndicator.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Options for adaptive indicators for p-adaptive DG scheme.
  \details   Options for adaptive indicators for p-adaptive DG scheme.
*/
// *****************************************************************************
#ifndef InciterPrefIndicatorOptions_h
#define InciterPrefIndicatorOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
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
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::pref_spectral_decay
                                  , kw::pref_non_conformity >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PrefIndicator() :
      tk::Toggle< PrefIndicatorType >(
        //! Group, i.e., options, name
        kw::pref_indicator::name(),
        //! Enums -> names
        { { PrefIndicatorType::SPECTRAL_DECAY,
            kw::pref_spectral_decay::name() },
          { PrefIndicatorType::NON_CONFORMITY,
            kw::pref_non_conformity::name() } },
        //! keywords -> Enums
        { { kw::pref_spectral_decay::string(),
            PrefIndicatorType::SPECTRAL_DECAY },
          { kw::pref_non_conformity::string(),
            PrefIndicatorType::NON_CONFORMITY } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterPrefIndicatorOptions_h
