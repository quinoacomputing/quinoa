//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Mon 29 Jul 2013 09:20:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that driver Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <Driver.h>

namespace Quinoa {

class Memory;
class Timer;
class Paradigm;
class Geometry;
class Physics;

//! QuinoaDriver : Driver
class QuinoaDriver : public Driver {

  public:
    //! Constructor
    QuinoaDriver(int argc,
                 char** argv,
                 Memory* const memory,
                 Paradigm* const paradigm);

    //! Destructor
    virtual ~QuinoaDriver() noexcept;

    //! Solve
    virtual void execute() const;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

  private:
    //! Don't permit copy constructor
    QuinoaDriver(const QuinoaDriver&) = delete;
    //! Don't permit assigment constructor
    QuinoaDriver& operator=(const QuinoaDriver&) = delete;
    //! Don't permit move constructor
    QuinoaDriver(QuinoaDriver&&) = delete;
    //! Don't permit move assignment
    QuinoaDriver& operator=(QuinoaDriver&&) = delete;

    //! Instantiate geometry object
    void initGeometry();

    //! Instantiate physics object
    void initPhysics();

    Memory* const m_memory;           //!< Memory object
    Paradigm* const m_paradigm;       //!< Parallel paradigm object

    Geometry* m_geometry;             //!< Geometry object
    Physics* m_physics;               //!< Physics object
    Timer* m_timer;                   //!< Timer object
};

} // namespace Quinoa

#endif // QuinoaDriver_h
