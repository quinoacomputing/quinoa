//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Sun 01 Jun 2014 11:45:27 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <make_unique.h>

#include <Factory.h>
#include <QuinoaDriver.h>
#include <Quinoa/InputDeck/Parser.h>
#include <Quinoa/CmdLine/Parser.h>
#include <AnalyticGeometry.h>
#include <DiscreteGeometry.h>
#include <HomMix.h>
#include <HomHydro.h>
#include <HomRT.h>
#include <SPINSFlow.h>
#include <TestSDE.h>

using quinoa::QuinoaDriver;

QuinoaDriver::QuinoaDriver(int argc, char** argv, const tk::Print& print)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \details   Instantiate geometry, MonteCarlo, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line into cmdline
  std::unique_ptr< ctr::CmdLine > cmdline;
  CmdLineParser cmdParser(argc, argv, print, cmdline);

  // Parse input deck into m_control, transfer cmdline (no longer needed)
  InputDeckParser inputdeckParser(print, std::move(cmdline), m_control);

  // Create essentials
  m_print = tk::make_unique< QuinoaPrint >( m_control );
  m_paradigm = tk::make_unique< tk::Paradigm >( print );
  m_timer = tk::make_unique< tk::Timer >();

  print.endpart();
  print.part("Factory");

  // Register random number generators
  initFactory( m_RNGFactory, m_paradigm->ompNthreads(),
               #ifdef HAS_MKL
               m_control->get< tag::param, tk::tag::rngmkl >(),
               #endif
               m_control->get< tag::param, tk::tag::rngsse >() );
  print.list< tk::ctr::RNG >
            ( "Registered random number generators", m_RNGFactory );

  // Bundle up essentials
  m_base = tk::make_unique< Base >( *m_print, *m_paradigm, *m_control,
                                    *m_timer, m_RNGFactory );

  //! Initialize factories
  initFactories( print );

  // Instantiate geometry
  ctr::GeometryType g = m_control->get< tag::selected, tag::geometry >();
  if (g != ctr::GeometryType::NO_GEOMETRY) {
    m_geometry = std::unique_ptr< Geometry >( m_geometryFactory[g]() );
  }

  // Instantiate MonteCarlo
  ctr::MonteCarloType m = m_control->get< tag::selected, tag::montecarlo >();
  if (m != ctr::MonteCarloType::NO_MONTECARLO) {
    m_montecarlo = std::unique_ptr< MonteCarlo >( m_MonteCarloFactory[m]() );
  }

  // Echo 'unspecified' if both geometry and MonteCarlo are unspecified
  if (!m_geometry && !m_montecarlo) {
    print.section( "Geometry", std::string("unspecified") );
    print.section( "MonteCarlo", std::string("unspecified") );
  }
}

void
QuinoaDriver::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize geometry and MonteCarlo factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register geometry types
  tk::record< AnalyticGeometry >
            ( m_geometryFactory, ctr::GeometryType::ANALYTIC, *m_base );
  tk::record< DiscreteGeometry >
            ( m_geometryFactory, ctr::GeometryType::DISCRETE, *m_base );
  print.list< ctr::Geometry >("Registered geometry types", m_geometryFactory);

  // Register MonteCarlo types
  tk::record< HomMix >
    ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_MIX, *m_base );
  tk::record< HomHydro >
    ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_HYDRO, *m_base );
  tk::record< HomRT >
    ( m_MonteCarloFactory, ctr::MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR,
      *m_base );
  tk::record< SPINSFlow >
    ( m_MonteCarloFactory, ctr::MonteCarloType::SPINSFLOW, *m_base );
  tk::record< TestSDE >
    ( m_MonteCarloFactory, ctr::MonteCarloType::TESTSDE, *m_base );
  print.list< ctr::MonteCarlo >
            ( "Registered MonteCarlo types", m_MonteCarloFactory );

  print.list( "Data layout policy",
              std::list< std::string > { ParProps(0,0).major(),
                                         "(Change CMake variable LAYOUT to "
                                         "change the data layout policy)" } );
}

void
QuinoaDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  //! Initialize and execute geometry (if any)
  if (m_geometry) {
    m_geometry->init();
    m_geometry->fill();
  }

  //! Run MonteCarlo (if any)
  if (m_montecarlo) {
    m_montecarlo->run();
  }
}
