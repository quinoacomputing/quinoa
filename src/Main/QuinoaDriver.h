//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Thu 07 Aug 2014 03:30:54 PM MDT
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
// 
//     //! Pointers to selected options
//     std::unique_ptr< MonteCarlo > m_montecarlo;  //!< MonteCarlo object

    const QuinoaPrint& m_print;
};

} // quinoa::

#endif // QuinoaDriver_h
