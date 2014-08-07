//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Thu 07 Aug 2014 03:30:42 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <make_unique.h>

#include <QuinoaDriver.h>
#include <Quinoa/InputDeck/Parser.h>
#include <quinoa.decl.h>
#include <Handler.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;

} // quinoa::

using quinoa::QuinoaDriver;

QuinoaDriver::QuinoaDriver( const QuinoaPrint& print,
                            const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \details   Instantiate MonteCarlo, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)
  try {

    // Parse input deck into g_inputdeck, transfer cmdline (no longer needed)
    InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );

    m_print.endpart();

  } catch (...) { tk::processException(); }

//   // Parse command line into cmdline
//   std::unique_ptr< ctr::CmdLine > cmdline;
//   CmdLineParser cmdParser(argc, argv, print, cmdline);
// 
//   // Parse input deck into m_control, transfer cmdline (no longer needed)
//   InputDeckParser inputdeckParser(print, std::move(cmdline), m_control);
// 
//   // Create essentials
//   m_print = tk::make_unique< QuinoaPrint >( m_control );
//   m_paradigm = tk::make_unique< tk::Paradigm >( print );
//   m_timer = tk::make_unique< tk::Timer >();
// 
//   print.endpart();
//   print.part("Factory");
// 
//   // Register random number generators
//   initFactory( m_RNGFactory, m_paradigm->ompNthreads(),
//                #ifdef HAS_MKL
//                m_control->get< tag::param, tk::tag::rngmkl >(),
//                #endif
//                m_control->get< tag::param, tk::tag::rngsse >() );
//   print.list< tk::ctr::RNG >
//             ( "Registered random number generators", m_RNGFactory );
// 
//   // Bundle up essentials
//   m_base = tk::make_unique< Base >( *m_print, *m_paradigm, *m_control,
//                                     *m_timer, m_RNGFactory );
// 
//   //! Initialize factories
//   initFactories( print );
// 
//   // Instantiate MonteCarlo
//   ctr::MonteCarloType m = m_control->get< tag::selected, tag::montecarlo >();
//   if (m != ctr::MonteCarloType::NO_MONTECARLO) {
//     m_montecarlo = std::unique_ptr< MonteCarlo >( m_MonteCarloFactory[m]() );
//   }
// 
//   // Echo 'unspecified' if MonteCarlo is unspecified
//   if (!m_montecarlo) {
//     print.section( "MonteCarlo", std::string("unspecified") );
//   }
}

void
QuinoaDriver::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize MonteCarlo factory
//! \author  J. Bakosi
//******************************************************************************
{
//   // Register MonteCarlo types
//   tk::record< HomMix >
//     ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_MIX, *m_base );
//   tk::record< HomHydro >
//     ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_HYDRO, *m_base );
//   tk::record< HomRT >
//     ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR,
//       *m_base );
//   tk::record< SPINSFlow >
//     ( m_MonteCarloFactory, ctr::MonteCarloType::SPINSFLOW, *m_base );
//   tk::record< TestSDE >
//     ( m_MonteCarloFactory, ctr::MonteCarloType::TESTSDE, *m_base );
//   print.list< ctr::MonteCarlo >
//             ( "Registered MonteCarlo types", m_MonteCarloFactory );
// 
//   print.list( "Data layout policy",
//               std::list< std::string > { ParProps(0,0).major(),
//                                          "(Change CMake variable LAYOUT to "
//                                          "change the data layout policy)" } );
}

void
QuinoaDriver::execute()
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
//   //! Run MonteCarlo (if any)
//   if (m_montecarlo) {
//     m_montecarlo->run();
//   }
}
