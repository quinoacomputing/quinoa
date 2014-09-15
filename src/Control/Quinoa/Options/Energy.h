//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Energy.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 08:13:11 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Energy model options and associations
  \details   Energy model options and associations
*/
//******************************************************************************
#ifndef QuinoaEnergyOptions_h
#define QuinoaEnergyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Energy model types
enum class EnergyType : uint8_t { NO_ENERGY=0 };

//! Pack/Unpack BatteryType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, EnergyType& e ) { tk::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class Energy : public tk::Toggle< EnergyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Energy() :
      Toggle< EnergyType >( "Energy",
      //! Enums -> names
      { { EnergyType::NO_ENERGY, "n/a" } },
      //! keywords -> Enums
      { { "no_energy", EnergyType::NO_ENERGY } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaEnergyOptions_h
