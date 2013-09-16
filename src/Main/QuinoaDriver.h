//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 10:19:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <Driver.h>
#include <Base.h>
#include <QuinoaControl.h>
#include <QuinoaPrint.h>

namespace quinoa {

class Geometry;
class Physics;

//! QuinoaDriver : Driver
class QuinoaDriver : public Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc,
                          char** argv,
                          Base& base);

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

    Base& m_base;                     //!< Essentials

    Geometry* m_geometry = nullptr;   //!< Geometry object
    Physics* m_physics = nullptr;     //!< Physics object
};

} // namespace quinoa

#endif // QuinoaDriver_h
