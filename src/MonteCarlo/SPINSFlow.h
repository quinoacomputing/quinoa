//******************************************************************************
/*!
  \file      src/MonteCarlo/SPINSFlow.h
  \author    J. Bakosi
  \date      Tue 31 Dec 2013 01:35:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************
#ifndef SPINSFlow_h
#define SPINSFlow_h

#include <Physics.h>

namespace quinoa {

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    explicit SPINSFlow(const Base& base) : Physics(base) {}

    //! Destructor
    ~SPINSFlow() override = default;

    //! Initialize
    void init() override;

    //! Run
    void run() override;

  private:
    //! Don't permit copy constructor
    SPINSFlow(const SPINSFlow&) = delete;
    //! Don't permit copy assigment
    SPINSFlow& operator=(const SPINSFlow&) = delete;
    //! Don't permit move constructor
    SPINSFlow(SPINSFlow&&) = delete;
    //! Don't permit move assigment
    SPINSFlow& operator=(SPINSFlow&&) = delete;
};

} // quinoa::

#endif // SPINSFlow_h
