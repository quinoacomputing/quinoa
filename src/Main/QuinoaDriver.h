//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Fri Oct 18 11:27:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <Driver.h>
#include <Base.h>
#include <Geometry.h>
#include <Physics.h>

//! Everything that contributes to the quina executable
namespace quinoa {

//! QuinoaDriver : Driver
class QuinoaDriver : public tk::Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc, char** argv, const tk::Print& print);

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

    //! Initialize geometry and physics factories
    void initFactory();

    std::unique_ptr< ctr::InputDeck > m_control; //!< Control
    std::unique_ptr< QuinoaPrint > m_print;      //!< Pretty printer
    std::unique_ptr< tk::Paradigm > m_paradigm;  //!< Parallel compute env.
    std::unique_ptr< tk::Timer > m_timer;        //!< Timer
    std::unique_ptr< Base > m_base;              //!< Essentials bundle

    //! Geometry and Physics factories
    std::map< ctr::PhysicsType, std::function<Physics*()> > m_physicsFactory;
    std::map< ctr::GeometryType, std::function<Geometry*()> > m_geometryFactory;

    std::unique_ptr< Geometry > m_geometry;   //!< Geometry object
    std::unique_ptr< Physics > m_physics;     //!< Physics object
};

} // quinoa::

#endif // QuinoaDriver_h
