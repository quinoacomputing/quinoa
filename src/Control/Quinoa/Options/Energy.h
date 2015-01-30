//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Energy.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 07:54:13 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Energy model options
  \details   Energy model options
*/
//******************************************************************************
#ifndef QuinoaEnergyOptions_h
#define QuinoaEnergyOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace quinoa {
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
} // quinoa::

#endif // QuinoaEnergyOptions_h
