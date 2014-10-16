//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Frequency.h
  \author    J. Bakosi
  \date      Tue 16 Sep 2014 08:16:44 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef QuinoaFrequencyOptions_h
#define QuinoaFrequencyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Frequency model types
enum class FrequencyType : uint8_t { NO_FREQUENCY=0,
                                     GAMMA };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, FrequencyType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class Frequency : public tk::Toggle< FrequencyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Frequency() :
      Toggle< FrequencyType >( "Turbulence frequency",
        //! Enums -> names
        { { FrequencyType::NO_FREQUENCY, "n/a" },
          { FrequencyType::GAMMA, kw::freq_gamma().name() } },
        //! keywords -> Enums
        { { "no_frequency", FrequencyType::NO_FREQUENCY },
          { kw::freq_gamma().string(), FrequencyType::GAMMA } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaFrequencyOptions_h
