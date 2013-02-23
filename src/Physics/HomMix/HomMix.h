//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Sat 23 Feb 2013 08:37:51 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mix model
  \details   Homogeneous material mix model
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

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    HomMix(Memory* const memory,
           Paradigm* const paradigm,
           Control* const control,
           Timer* const timer);

    //! Destructor
    virtual ~HomMix();

    //! Echo informaion on model
    virtual void echo();

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
    void outJPDF(const real t);

    //! One-liner report
    void reportHeader();
    void report(const int it,
                const int nstep,
                const real t,
                const real dt,
                const bool wroteJPDF);

    const int m_nscalar;                //!< Number of mixing scalars
    const real m_term;                  //!< Maximum time to simulate
    const string m_jpdf_filename_base;  //!< Joint PDF filename base
    const TimerIdx m_totalTime;         //!< Timer measuring the total run
    Mix* m_mix;                         //!< Mix model object
};

} // namespace Quinoa

#endif // HomMix_h
