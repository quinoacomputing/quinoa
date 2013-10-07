//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:38:23 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <Physics.h>
#include <Base.h>

namespace quinoa {

class Memory;
class Paradigm;
class QuinoaControl;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    explicit HomMix(const Base& base);

    //! Destructor
    ~HomMix() override = default;

    //! Initialize model
    void init() override;

    //! Solve model
    void solve() override;

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

    //! Echo information on homogeneous material mixing physics
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

#endif // HomMix_h
