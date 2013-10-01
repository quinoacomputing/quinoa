//******************************************************************************
/*!
  \file      src/Physics/HomHydro/HomHydro.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 08:39:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous hydrodynamics
  \details   Homogeneous hydrodynamics
*/
//******************************************************************************
#ifndef HomHydro_h
#define HomHydro_h

#include <Physics.h>
#include <Hydro/Hydro.h>
#include <Timer.h>

namespace quinoa {

class Memory;
class Paradigm;
class Statistics;
class GlobWriter;
class TxtStatWriter;

//! HomHydro : Physics
class HomHydro : public Physics {

  public:
    //! Constructor
    explicit HomHydro(const Base& base);

    //! Destructor
    ~HomHydro() override = default;

    //! Initialize model
    void init() override;

    //! Solve model
    void solve() override;

  private:
    //! Don't permit copy constructor
    HomHydro(const HomHydro&) = delete;
    //! Don't permit copy assigment
    HomHydro& operator=(const HomHydro&) = delete;
    //! Don't permit move constructor
    HomHydro(HomHydro&&) = delete;
    //! Don't permit move assigment
    HomHydro& operator=(HomHydro&&) = delete;
    
    //! Echo information on homogeneous hydrodynamics physics
    void echo();

    //! One-liner report
    void reportHeader();
    void report(const uint64_t it,
                const uint64_t nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wroteStat);

    //! Advance
    void advance(real dt);

    //! Output joint PDF
    void outJpdf(const real t);

    const TimerIdx m_totalTime;           //!< Timer measuring the total run
};

} // namespace quinoa

#endif // HomHydro_h
