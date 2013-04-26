//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Fri Apr 26 17:06:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <map>

#include <Physics.h>
#include <ControlTypes.h>
#include <Timer.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Mix;
class Timer;
class Statistics;
class GlobWriter;
class TxtPlotWriter;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    explicit HomMix(Memory* const memory,
                    Paradigm* const paradigm,
                    Control* const control,
                    Timer* const timer);

    //! Destructor
    virtual ~HomMix() noexcept;

    //! Echo informaion on model
    virtual void echo() const;

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

    //! Output joint scalar PDF
    void outJpdf(const real t);

    //! One-liner report
    void reportHeader() const;
    void report(const int it,
                const int nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    const int m_nscalar;                  //!< Number of mixing scalars
    const real m_term;                    //!< Maximum time to simulate
    const TimerIdx m_totalTime;           //!< Timer measuring the total run

    Mix* m_mix;                           //!< Mix model object
    Statistics* m_statistics;             //!< Statistics estimator object
    GlobWriter* m_glob;                   //!< Glob file writer
    TxtPlotWriter* m_plot;                //!< Plot file writer
};

} // namespace Quinoa

#endif // HomMix_h
