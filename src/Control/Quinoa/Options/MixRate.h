//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MixRate.h
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 08:00:24 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Turbulence frequency model options
  \details   Turbulence frequency model options
*/
//******************************************************************************
#ifndef QuinoaMixRateOptions_h
#define QuinoaMixRateOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace quinoa {
namespace ctr {

//! Material mix rate model types
//! \author J. Bakosi
enum class MixRateType : uint8_t { NO_MIXRATE=0,
                                   GAMMA };

//! \brief Pack/Unpack MixRateType: forward overload to generic enum class
//!   packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, MixRateType& e ) { PUP::pup( p, e ); }

//! \brief Mix rate model options: outsource searches to base templated on
//!   enum type
//! \author J. Bakosi
class MixRate : public tk::Toggle<MixRateType> {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::mixrate_gamma
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit MixRate() :
      Toggle< MixRateType >(
        //! Group, i.e., options, name
        "Material mix rate",
        //! Enums -> names
        { { MixRateType::NO_MIXRATE, "n/a" },
          { MixRateType::GAMMA, kw::mixrate_gamma::name() } },
        //! keywords -> Enums
        { { "no_mixrate", MixRateType::NO_MIXRATE },
          { kw::mixrate_gamma::string(), MixRateType::GAMMA } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaMixRateOptions_h
