//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:43:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************
#ifndef SPINSFlow_h
#define SPINSFlow_h

#include <Physics.h>
#include <Hydro/Hydro.h>
#include <Timer.h>

namespace quinoa {

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    explicit SPINSFlow(const Base& base) :
      Physics(base),
      m_totalTime(base.timer.create("Total solution")) {}

    //! Destructor
    ~SPINSFlow() override = default;

    //! Initialize model
    void init() override;

    //! Solve model
    void solve() override;

  private:
    //! Don't permit copy constructor
    SPINSFlow(const SPINSFlow&) = delete;
    //! Don't permit copy assigment
    SPINSFlow& operator=(const SPINSFlow&) = delete;
    //! Don't permit move constructor
    SPINSFlow(SPINSFlow&&) = delete;
    //! Don't permit move assigment
    SPINSFlow& operator=(SPINSFlow&&) = delete;

    //! Echo information on standalone-particle incompressible Navier-Stokes
    //! physics
    void echo();

    const tk::TimerIdx m_totalTime;           //!< Timer measuring the total run
};

} // quinoa::

#endif // SPINSFlow_h
