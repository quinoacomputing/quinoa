//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MonteCarlo.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 11:04:07 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     MonteCarlo options and associations
  \details   MonteCarlo options and associations
*/
//******************************************************************************
#ifndef QuinoaMonteCarloOptions_h
#define QuinoaMonteCarloOptions_h

#include <map>

#include <Toggle.h>
#include <Keywords.h>

namespace quinoa {
namespace ctr {

//! MonteCarlo types
enum class MonteCarloType : uint8_t { NO_MONTECARLO=0,
                                      HOMOGENEOUS_MIX,
                                      HOMOGENEOUS_HYDRO,
                                      HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                      SPINSFLOW };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, MonteCarloType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class MonteCarlo : public tk::Toggle< MonteCarloType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MonteCarlo() :
      Toggle< MonteCarloType >( "MonteCarlo",
        //! Enums -> names
        { { MonteCarloType::NO_MONTECARLO, "n/a" },
          { MonteCarloType::HOMOGENEOUS_MIX, kw::hommix().name() },
          { MonteCarloType::HOMOGENEOUS_HYDRO, kw::homhydro().name() },
          { MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR, kw::homrt().name() },
          { MonteCarloType::SPINSFLOW, kw::spinsflow().name() } },
        //! keywords -> Enums
        { { "no_montecarlo", MonteCarloType::NO_MONTECARLO },
          { kw::hommix().string(), MonteCarloType::HOMOGENEOUS_MIX },
          { kw::homhydro().string(), MonteCarloType::HOMOGENEOUS_HYDRO },
          { kw::homrt().string(), MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR },
          { kw::spinsflow().string(), MonteCarloType::SPINSFLOW } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaMonteCarloOptions_h
