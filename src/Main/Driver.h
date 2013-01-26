//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 08:04:29 AM MST
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
class Control;

//! Driver base class
class Driver {

  public:
    //! Constructor
    Driver(int argc, char** argv, Memory* memory, Paradigm* paradigm);

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

    Memory* m_memory;                 //!< Memory object
    Paradigm* m_paradigm;             //!< Parallel paradigm object
    Physics* m_physics;               //!< Physics object
    Control* m_control;               //!< Control object

    int m_argc;                       //!< Argument count from command line
    char** m_argv;                    //!< Argument vector from command line
};

} // namespace Quinoa

#endif // Driver_h
