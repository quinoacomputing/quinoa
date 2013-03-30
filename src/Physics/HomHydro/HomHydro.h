//******************************************************************************
/*!
  \file      src/Physics/HomHydro/HomHydro.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 11:43:14 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous hydrodynamics
  \details   Homogeneous hydrodynamics
*/
//******************************************************************************
#ifndef HomHydro_h
#define HomHydro_h

#include <map>

#include <Physics.h>
#include <ControlTypes.h>
#include <Timer.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Hydro;
class Timer;
class Statistics;
class GlobWriter;
class TxtPlotWriter;

//! HomHydro : Physics
class HomHydro : public Physics {

  public:
    //! Constructor
    HomHydro(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
             Timer* const timer);

    //! Destructor
    virtual ~HomHydro();

    //! Echo informaion on model
    virtual void echo();

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

    //! Output joint PDF
    void outJpdf(const real t);

    //! One-liner report
    void reportHeader();
    void report(const int it,
                const int nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    const real m_term;                    //!< Maximum time to simulate
    const TimerIdx m_totalTime;           //!< Timer measuring the total run

    Hydro* m_hydro;                       //!< Hydro model object
    Statistics* m_statistics;             //!< Statistics estimator object
    GlobWriter* m_glob;                   //!< Glob file writer
    TxtPlotWriter* m_plot;                //!< Plot file writer
};

} // namespace Quinoa

#endif // HomHydro_h
