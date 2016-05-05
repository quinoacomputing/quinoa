//******************************************************************************
/*!
  \file      src/Control/Breeze/Options/MonteCarlo.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:21:41 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     MonteCarlo options
  \details   MonteCarlo options
*/
//******************************************************************************
#ifndef BreezeMonteCarloOptions_h
#define BreezeMonteCarloOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace breeze {
namespace ctr {

//! MonteCarlo types
//! \author J. Bakosi
enum class MonteCarloType : uint8_t { NO_MONTECARLO=0,
                                      HOMOGENEOUS_MIX,
                                      HOMOGENEOUS_HYDRO,
                                      HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                      SPINSFLOW };

//! \brief Pack/Unpack MonteCarloType: forward overload to generic enum class
//!    packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, MonteCarloType& e ) { PUP::pup( p, e ); }

//! \brief MonteCarlo options: outsource searches to base templated on enum type
//! \author J. Bakosi
class MonteCarlo : public tk::Toggle< MonteCarloType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::hommix
                                       , kw::homhydro
                                       , kw::homrt
                                       , kw::spinsflow
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit MonteCarlo() :
      Toggle< MonteCarloType >(
        //! Group, i.e., options, name
        "MonteCarlo",
        //! Enums -> names
        { { MonteCarloType::NO_MONTECARLO, "n/a" },
          { MonteCarloType::HOMOGENEOUS_MIX, kw::hommix::name() },
          { MonteCarloType::HOMOGENEOUS_HYDRO, kw::homhydro::name() },
          { MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR, kw::homrt::name() },
          { MonteCarloType::SPINSFLOW, kw::spinsflow::name() } },
        //! keywords -> Enums
        { { "no_montecarlo", MonteCarloType::NO_MONTECARLO },
          { kw::hommix::string(), MonteCarloType::HOMOGENEOUS_MIX },
          { kw::homhydro::string(), MonteCarloType::HOMOGENEOUS_HYDRO },
          { kw::homrt::string(), MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR },
          { kw::spinsflow::string(), MonteCarloType::SPINSFLOW } } ) {}
};

} // ctr::
} // breeze::

#endif // BreezeMonteCarloOptions_h
