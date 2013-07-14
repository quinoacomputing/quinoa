//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Sat 13 Jul 2013 08:32:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class declaration
  \details   Driver base class declaration
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <Memory.h>

namespace Quinoa {

class Geometry;
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
    ~Driver() noexcept;

    //! Solve
    void execute() const;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;

    //! Instantiate geometry object
    void initGeometry();

    //! Instantiate physics object
    void initPhysics();

    Memory* const m_memory;           //!< Memory object
    Paradigm* const m_paradigm;       //!< Parallel paradigm object

    Geometry* m_geometry;             //!< Geometry object
    Physics* m_physics;               //!< Physics object
    Control* m_control;               //!< Control object
    Timer* m_timer;                   //!< Timer object
};

} // namespace Quinoa

#endif // Driver_h
