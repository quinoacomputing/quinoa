//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:08:07 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <QuinoaPrint.h>
#include <Quinoa/CmdLine/CmdLine.h>

//! Everything that contributes to the quinoa executable
namespace quinoa {

// //! Geometry factory type
// using GeometryFactory =
//   std::map< ctr::GeometryType, std::function< Geometry*() > >;
// 
// //! MonteCarlo factory type
// using MonteCarloFactory =
//   std::map< ctr::MonteCarloType, std::function< MonteCarlo*() > >;

//! Quinoa driver used polymorphically with Driver
class QuinoaDriver {

  public:
    //! Constructor
    explicit QuinoaDriver( const QuinoaPrint& print,
                           const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute();

  private:
    //! Initialize factories
    void initFactories( const tk::Print& print );

//     //! Factories
//     tk::RNGFactory m_RNGFactory;                 //!< RNG factory
//     MonteCarloFactory m_MonteCarloFactory;       //!< MonteCarlo factory
//     GeometryFactory m_geometryFactory;           //!< Geometry factory
// 
//     //! Pointers to selected options
//     std::unique_ptr< Geometry > m_geometry;      //!< Geometry object
//     std::unique_ptr< MonteCarlo > m_montecarlo;  //!< MonteCarlo object

    const QuinoaPrint& m_print;
};

} // quinoa::

#endif // QuinoaDriver_h
