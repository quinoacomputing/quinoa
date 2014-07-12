//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MonteCarlo.h
  \author    J. Bakosi
  \date      Wed 29 Jan 2014 09:34:55 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MonteCarlo options and associations
  \details   MonteCarlo options and associations
*/
//******************************************************************************
#ifndef QuinoaMonteCarloOptions_h
#define QuinoaMonteCarloOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! MonteCarlo types
enum class MonteCarloType : uint8_t { NO_MONTECARLO=0,
                                      TESTSDE,
                                      HOMOGENEOUS_MIX,
                                      HOMOGENEOUS_HYDRO,
                                      HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                      SPINSFLOW };

//! Class with base templated on the above enum class with associations
class MonteCarlo : public tk::Toggle<MonteCarloType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MonteCarlo() : Toggle<MonteCarloType>("MonteCarlo", names, values)
    {}

  private:
    //! Don't permit copy constructor
    MonteCarlo(const MonteCarlo&) = delete;
    //! Don't permit copy assigment
    MonteCarlo& operator=(const MonteCarlo&) = delete;
    //! Don't permit move constructor
    MonteCarlo(MonteCarlo&&) = delete;
    //! Don't permit move assigment
    MonteCarlo& operator=(MonteCarlo&&) = delete;

    //! Get access to Monte Carlo keywords
    const kw::testsde testsde {};
    const kw::hommix hommix {};
    const kw::homhydro homhydro {};
    const kw::homrt homrt {};
    const kw::spinsflow spinsflow {};

    //! Enums -> names
    const std::map<MonteCarloType, std::string> names {
      { MonteCarloType::NO_MONTECARLO, "n/a" },
      { MonteCarloType::TESTSDE, testsde.name() },
      { MonteCarloType::HOMOGENEOUS_MIX, hommix.name() },
      { MonteCarloType::HOMOGENEOUS_HYDRO, homhydro.name() },
      { MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR, homrt.name() },
      { MonteCarloType::SPINSFLOW, spinsflow.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, MonteCarloType> values {
      { "no_montecarlo", MonteCarloType::NO_MONTECARLO },
      { testsde.string(), MonteCarloType::TESTSDE },
      { hommix.string(), MonteCarloType::HOMOGENEOUS_MIX },
      { homhydro.string(), MonteCarloType::HOMOGENEOUS_HYDRO },
      { homrt.string(), MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR },
      { spinsflow.string(), MonteCarloType::SPINSFLOW }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaMonteCarloOptions_h
