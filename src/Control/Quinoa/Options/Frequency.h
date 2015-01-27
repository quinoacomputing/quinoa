//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Frequency.h
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 08:00:05 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Turbulence frequency model options
  \details   Turbulence frequency model options
*/
//******************************************************************************
#ifndef QuinoaFrequencyOptions_h
#define QuinoaFrequencyOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace quinoa {
namespace ctr {

//! Frequency model types
//! \author J. Bakosi
enum class FrequencyType : uint8_t { NO_FREQUENCY=0,
                                     GAMMA };

//! \brief Pack/Unpack FrequencyType: forward overload to generic enum class
//!    packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, FrequencyType& e ) { PUP::pup( p, e ); }

//! \brief Frequency model options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class Frequency : public tk::Toggle< FrequencyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::freq_gamma
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Frequency() :
      Toggle< FrequencyType >(
        //! Group, i.e., options, name
        "Turbulence frequency",
        //! Enums -> names
        { { FrequencyType::NO_FREQUENCY, "n/a" },
          { FrequencyType::GAMMA, kw::freq_gamma::name() } },
        //! keywords -> Enums
        { { "no_frequency", FrequencyType::NO_FREQUENCY },
          { kw::freq_gamma::string(), FrequencyType::GAMMA } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaFrequencyOptions_h
