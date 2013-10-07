//******************************************************************************
/*!
  \file      src/Physics/HomRT/HomRT.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:40:47 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Rayleigh-Taylor
  \details   Homogeneous Rayleigh-Taylor
*/
//******************************************************************************
#ifndef HomRT_h
#define HomRT_h

#include <map>

#include <Physics.h>
#include <Mix/Mix.h>
#include <Timer.h>

namespace quinoa {

//! HomRT : Physics
class HomRT : public Physics {

  public:
    //! Constructor
    explicit HomRT(const Base& base);

    //! Destructor
    ~HomRT() override = default;

    //! Initialize model
    void init() override;

    //! Solve model
    void solve() override;

  private:
    //! Don't permit copy constructor
    HomRT(const HomRT&) = delete;
    //! Don't permit copy assigment
    HomRT& operator=(const HomRT&) = delete;
    //! Don't permit move constructor
    HomRT(HomRT&&) = delete;
    //! Don't permit move assigment
    HomRT& operator=(HomRT&&) = delete;

    //! Echo information on homogeneous Rayleigh-Taylor physics
    void echo();

    //! One-liner report
    void reportHeader() const;
    void report(const uint64_t it,
                const uint64_t nstep,
                const tk::real t,
                const tk::real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    //! Advance
    void advance(tk::real dt);

    //! Output joint scalar PDF
    void outJpdf(const tk::real t);

    const tk::TimerIdx m_totalTime;           //!< Timer measuring the total run    
};

} // quinoa::

#endif // HomRT_h
