//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Wed Aug 28 15:10:16 2013
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
class Paradigm;
class Geometry;
class Control;
class Physics;

//! QuinoaDriver : Driver
class QuinoaDriver : public Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc,
                          char** argv,
                          Memory* const memory,
                          Paradigm* const paradigm);

    //! Destructor
    virtual ~QuinoaDriver() noexcept;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    virtual void finalize() noexcept;

    //! Solve
    virtual void execute() const;

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

    Control* m_control;               //!< Control object
    Geometry* m_geometry;             //!< Geometry object
    Physics* m_physics;               //!< Physics object
};

} // namespace Quinoa

#endif // QuinoaDriver_h
