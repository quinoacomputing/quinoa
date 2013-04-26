//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Fri Apr 26 12:43:23 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class declaration
  \details   Driver base class declaration
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <Memory.h>

namespace Quinoa {

class Physics;
class Timer;
class Control;

//! Driver base class
class Driver {

  public:
    //! Constructor
    Driver(int argc,
           char** argv,
           Memory* const memory,
           Paradigm* const paradigm);

    //! Destructor
    ~Driver();

    //! Setup: instantiate model, set initial conditions
    void setup();

    //! Solve
    void solve();

    //! Finalize (either at normal exit, or due to exception)
    void finalize();

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;

    Memory* const m_memory;           //!< Memory object
    Paradigm* const m_paradigm;       //!< Parallel paradigm object

    int m_argc;                       //!< Argument count from command line
    char** m_argv;                    //!< Argument vector from command line
    Physics* m_physics;               //!< Physics object
    Control* m_control;               //!< Control object
    Timer* m_timer;                   //!< Timer object
};

} // namespace Quinoa

#endif // Driver_h
