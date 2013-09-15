//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Sat 14 Sep 2013 08:15:55 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <Driver.h>
#include <QuinoaControl.h>
#include <QuinoaPrinter.h>

namespace quinoa {

class Memory;
class Paradigm;
class Geometry;
class QuinoaControl;
class Physics;

//! QuinoaDriver : Driver
class QuinoaDriver : public Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc,
                          char** argv,
                          Memory* const memory,
                          Paradigm* const paradigm,
                          const QuinoaPrinter& print);

    //! Destructor
    ~QuinoaDriver() noexcept override;

    //! Solve
    void execute() const override;

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
    const QuinoaPrinter& m_print;     //!< Pretty printer object

    QuinoaControl m_control;          //!< Control object
    Geometry* m_geometry;             //!< Geometry object
    Physics* m_physics;               //!< Physics object
};

} // namespace quinoa

#endif // QuinoaDriver_h
