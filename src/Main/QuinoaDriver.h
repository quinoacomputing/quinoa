//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:24:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <Driver.h>
#include <Base.h>
#include <QuinoaPrint.h>
#include <Geometry.h>
#include <Physics.h>

namespace quinoa {

//! QuinoaDriver : Driver
class QuinoaDriver : public Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc, char** argv, Base& base);

    //! Destructor
    ~QuinoaDriver() noexcept override = default;

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

    //! Initialize physics factory
    void initFactory();

    Base& m_base;                           //!< Essentials

    //! Geometry and Physics factories
    std::map<sel::PhysicsType, std::function<Physics*()>> m_physicsFactory;
    std::map<sel::GeometryType, std::function<Geometry*()>> m_geometryFactory;

    std::unique_ptr<Geometry> m_geometry;   //!< Geometry object
    std::unique_ptr<Physics> m_physics;     //!< Physics object
};

} // namespace quinoa

#endif // QuinoaDriver_h
