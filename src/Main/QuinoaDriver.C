//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Sat 01 Feb 2014 10:53:29 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

  // Create pretty printer
  m_print = tk::make_unique< QuinoaPrint >( m_control );
  m_paradigm = tk::make_unique< tk::Paradigm >( print );
  m_timer = tk::make_unique< tk::Timer >();

  print.endpart();
  print.part("Factory");

  // Register random number generators
  tk::ctr::RNG rng;
  std::list< tk::ctr::RNGType > regRNG;
  initRNGFactory( m_RNGFactory, rng, regRNG, m_paradigm->nthreads(),
                  m_control->get< tag::param, tk::tag::mklrng >(),
                  m_control->get< tag::param, tk::tag::rngsse >() );
  print.list("Registered random number generators", rng, regRNG);

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
  // Register geometries
  ctr::Geometry geometry;
  std::list< ctr::GeometryType > regGeo;
  tk::regist< AnalyticGeometry >( m_geometryFactory, regGeo, geometry,
                                 ctr::GeometryType::ANALYTIC, *m_base );
  tk::regist< DiscreteGeometry >( m_geometryFactory, regGeo, geometry,
                                  ctr::GeometryType::DISCRETE, *m_base );
  print.list("Registered geometries", geometry, regGeo);

  // Register MonteCarlo
  ctr::MonteCarlo mc;
  std::list< ctr::MonteCarloType > regMC;
  tk::regist< HomMix >( m_MonteCarloFactory, regMC, mc,
                        ctr::MonteCarloType::HOMOGENEOUS_MIX,
                        *m_base );
  tk::regist< HomHydro >( m_MonteCarloFactory, regMC, mc,
                          ctr::MonteCarloType::HOMOGENEOUS_HYDRO,
                          *m_base );
  tk::regist< HomRT >( m_MonteCarloFactory, regMC, mc,
                       ctr::MonteCarloType::HOMOGENEOUS_RAYLEIGH_TAYLOR,
                       *m_base );
  tk::regist< SPINSFlow >( m_MonteCarloFactory, regMC, mc,
                           ctr::MonteCarloType::SPINSFLOW,
                           *m_base );
  tk::regist< TestSDE >( m_MonteCarloFactory, regMC, mc,
                           ctr::MonteCarloType::TESTSDE,
                           *m_base );
  print.list("Registered MonteCarlo", mc, regMC);

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
