//******************************************************************************
/*!
  \file      src/Control/Breeze/Options/Energy.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:21:08 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Energy model options
  \details   Energy model options
*/
//******************************************************************************
#ifndef BreezeEnergyOptions_h
#define BreezeEnergyOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace breeze {
namespace ctr {

//! Energy model types
//! \author J. Bakosi
enum class EnergyType : uint8_t { NO_ENERGY=0 };

//! Pack/Unpack EnergyType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, EnergyType& e ) { PUP::pup( p, e ); }

//! \brief Energy model options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class Energy : public tk::Toggle< EnergyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector<
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Energy() :
      Toggle< EnergyType >(
        //! Group, i.e., options, name
        "Energy",
        //! Enums -> names
        { { EnergyType::NO_ENERGY, "n/a" } },
        //! keywords -> Enums
        { { "no_energy", EnergyType::NO_ENERGY } } ) {}
};

} // ctr::
} // breeze::

#endif // BreezeEnergyOptions_h
