// *****************************************************************************
/*!
  \file      src/Control/Breeze/Options/Mass.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:20:59 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mass model options
  \details   Mass model options
*/
// *****************************************************************************
#ifndef BreezeMassOptions_h
#define BreezeMassOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace breeze {
namespace ctr {

//! Mass model types
//! \author J. Bakosi
enum class MassType : uint8_t { NO_MASS=0,
                                BETA };

//! \brief Pack/Unpack MassType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, MassType& e ) { PUP::pup( p, e ); }

//! \brief Mass model options: outsource searches to base templated on enum type
//! \author J. Bakosi
class Mass : public tk::Toggle< MassType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::beta
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Mass() :
      Toggle< MassType >(
        //! Group, i.e., options, name
        "Mass",
        //! Enums -> names
        { { MassType::NO_MASS, "n/a" },
          { MassType::BETA, kw::mass_beta::name() } },
        //! keywords -> Enums
        { { "no_mass", MassType::NO_MASS },
          { kw::mass_beta::string(), MassType::BETA } } ) {}
};

} // ctr::
} // breeze::

#endif // BreezeMassOptions_h
