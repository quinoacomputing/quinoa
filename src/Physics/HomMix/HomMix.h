//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 09:05:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <map>

#include <Physics.h>
#include <Timer.h>

namespace Quinoa {

class Memory;
class Paradigm;
class QuinoaControl;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    explicit HomMix(Memory* const memory,
                    Paradigm* const paradigm,
                    QuinoaControl* const control,
                    Timer* const timer);

    //! Destructor
    virtual ~HomMix() noexcept = default;

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

    //! One-liner report
    void reportHeader() const;
    void report(const uint64_t it,
                const uint64_t nstep,
                const real t,
                const real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    //! Advance
    void advance(real dt);

    //! Output joint scalar PDF
    void outJpdf(const real t);

    const TimerIdx m_totalTime;           //!< Timer measuring the total run    
};

} // namespace Quinoa

#endif // HomMix_h
