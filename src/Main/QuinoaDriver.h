//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Wed Mar 19 07:58:08 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <RNGDriver.h>
#include <Base.h>
#include <Geometry.h>
#include <Physics.h>

//! Everything that contributes to the quina executable
namespace quinoa {

//! Geometry factory type
using GeometryFactory =
  std::map< ctr::GeometryType, std::function< Geometry*() > >;

//! MonteCarlo factory type
using MonteCarloFactory =
  std::map< ctr::MonteCarloType, std::function< MonteCarlo*() > >;

//! QuinoaDriver : Driver
class QuinoaDriver : public tk::RNGDriver {

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
    MonteCarloFactory m_MonteCarloFactory;       //!< MonteCarlo factory
    GeometryFactory m_geometryFactory;           //!< Geometry factory

    //! Pointers to selected options
    std::unique_ptr< Geometry > m_geometry;      //!< Geometry object
    std::unique_ptr< MonteCarlo > m_montecarlo;  //!< MonteCarlo object

    //! Pointers to essentials to be created
    std::unique_ptr< ctr::InputDeck > m_control; //!< Control
    std::unique_ptr< QuinoaPrint > m_print;      //!< Pretty printer
    std::unique_ptr< tk::Paradigm > m_paradigm;  //!< Parallel compute env.
    std::unique_ptr< tk::Timer > m_timer;        //!< Timer
    std::unique_ptr< Base > m_base;              //!< Essentials bundle
};

} // quinoa::

#endif // QuinoaDriver_h
