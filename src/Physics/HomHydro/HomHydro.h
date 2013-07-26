//******************************************************************************
/*!
  \file      src/Physics/HomHydro/HomHydro.h
  \author    J. Bakosi
  \date      Fri Jul 26 12:49:02 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous hydrodynamics
  \details   Homogeneous hydrodynamics
*/
//******************************************************************************
#ifndef HomHydro_h
#define HomHydro_h

#include <Physics.h>
#include <Hydro.h>
#include <Timer.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Statistics;
class GlobWriter;
class TxtStatWriter;

//! HomHydro : Physics
class HomHydro : public Physics {

  public:
    //! Constructor
    explicit HomHydro(Memory* const memory,
                      Paradigm* const paradigm,
                      Control* const control,
                      Timer* const timer);

    //! Destructor
    virtual ~HomHydro() noexcept = default;

    //! Echo informaion on model
    virtual void echo() const;

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomHydro(const HomHydro&) = delete;
    //! Don't permit copy assigment
    HomHydro& operator=(const HomHydro&) = delete;
    //! Don't permit move constructor
    HomHydro(HomHydro&&) = delete;
    //! Don't permit move assigment
    HomHydro& operator=(HomHydro&&) = delete;

    //! One-liner report
    void reportHeader();
    void report(const uint64_t it,
                const uint64_t nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wroteStat);

    //! Output joint PDF
    void outJpdf(const real t);

    const TimerIdx m_totalTime;           //!< Timer measuring the total run
};

} // namespace Quinoa

#endif // HomHydro_h
