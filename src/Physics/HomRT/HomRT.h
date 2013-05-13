//******************************************************************************
/*!
  \file      src/Physics/HomRT/HomRT.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 08:15:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Rayleigh-Taylor
  \details   Homogeneous Rayleigh-Taylor
*/
//******************************************************************************
#ifndef HomRT_h
#define HomRT_h

#include <map>

#include <Physics.h>
#include <Mix.h>
#include <Timer.h>
#include <ControlTypes.h>

namespace Quinoa {

class Memory;
class Paradigm;

//! HomRT : Physics
class HomRT : public Physics {

  public:
    //! Constructor
    explicit HomRT(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   Timer* const timer);

    //! Destructor
    virtual ~HomRT() noexcept = default;

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomRT(const HomRT&) = delete;
    //! Don't permit copy assigment
    HomRT& operator=(const HomRT&) = delete;
    //! Don't permit move constructor
    HomRT(HomRT&&) = delete;
    //! Don't permit move assigment
    HomRT& operator=(HomRT&&) = delete;

    //! One-liner report
    void reportHeader() const;
    void report(const int it,
                const int nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    //! Output joint scalar PDF
    void outJpdf(const real t);

    const TimerIdx m_totalTime;           //!< Timer measuring the total run    
};

} // namespace Quinoa

#endif // HomRT_h
