// *****************************************************************************
/*!
  \file      src/Control/RNGTest/Options/Battery.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test suite batteries options
  \details   Random number generator test suite batteries options
*/
// *****************************************************************************
#ifndef RNGTestBatteryOptions_h
#define RNGTestBatteryOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace rngtest {
namespace ctr {

//! Random number generator battery types
//! \author J. Bakosi
enum class BatteryType : uint8_t { NO_BATTERY=0,
                                   SMALLCRUSH,
                                   CRUSH,
                                   BIGCRUSH };

//! Pack/Unpack BatteryType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, BatteryType& e ) { PUP::pup( p, e ); }

//! \brief Battery options: outsource searches to base templated on enum type
//! \author J. Bakosi
class Battery : public tk::Toggle< BatteryType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::smallcrush
                                       , kw::crush
                                       , kw::bigcrush
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Battery() :
      tk::Toggle< BatteryType >(
        //! Group, i.e., options, name
        "RNG battery",
        //! Enums -> names
        { { BatteryType::NO_BATTERY, "n/a" },
          { BatteryType::SMALLCRUSH, kw::smallcrush::name() },
          { BatteryType::CRUSH, kw::crush::name() },
          { BatteryType::BIGCRUSH, kw::bigcrush::name() } },
        //! keywords -> Enums
        { { "no_battery", BatteryType::NO_BATTERY },
          { kw::smallcrush::string(), BatteryType::SMALLCRUSH },
          { kw::crush::string(), BatteryType::CRUSH },
          { kw::bigcrush::string(), BatteryType::BIGCRUSH } } ) {}
};

} // ctr::
} // rngtest::

#endif // RNGTestBatteryOptions_h
