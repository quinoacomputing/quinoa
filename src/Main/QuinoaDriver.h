//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Thu Nov 14 09:48:42 2013
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

//! Geometry factory type
using GeometryFactory = std::map< ctr::GeometryType, std::function< Geometry*() > >;

//! Physics factory type
using PhysicsFactory = std::map< ctr::PhysicsType, std::function< Physics*() > >;

//! QuinoaDriver : Driver
class QuinoaDriver : public tk::Driver {

  public:
    //! Constructor
    explicit QuinoaDriver(int argc, char** argv, const tk::Print& print);

    //! Destructor
    ~QuinoaDriver() noexcept override = default;

    //! Execute
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

    //! Initialize factories
    void initFactories(const tk::Print& print);

    //! Factories
    tk::RNGFactory m_RNGFactory;                 //!< RNG factory
    PhysicsFactory m_physicsFactory;             //!< Physics factory
    GeometryFactory m_geometryFactory;           //!< Geometry factory

    //! Pointers to selected options
    std::unique_ptr< Geometry > m_geometry;      //!< Geometry object
    std::unique_ptr< Physics > m_physics;        //!< Physics object

    //! Pointers to essentials to be created
    std::unique_ptr< ctr::InputDeck > m_control; //!< Control
    std::unique_ptr< QuinoaPrint > m_print;      //!< Pretty printer
    std::unique_ptr< tk::Paradigm > m_paradigm;  //!< Parallel compute env.
    std::unique_ptr< tk::Timer > m_timer;        //!< Timer
    std::unique_ptr< Base > m_base;              //!< Essentials bundle
};

} // quinoa::

#endif // QuinoaDriver_h
