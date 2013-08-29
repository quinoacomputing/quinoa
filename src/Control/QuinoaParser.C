//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:31:21 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa control file parser
  \details   Quinoa control file parser
*/
//******************************************************************************

#include <pegtl.hh>

#include <QuinoaParser.h>
#include <QuinoaGrammar.h>
#include <Control.h>

using namespace quinoa;

void
QuinoaParser::parse()
//******************************************************************************
//  Parse quinoa control file
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize new bundle for parsed data with defaults
  control::Bundle stack(control::defaults);
  // Initialize new bool bundle for indicating what data is set in bundle
  control::BoolBundle
    boolstack(std::tuple_size<decltype(control::defaults)>::value, false);

  //std::cout << "==== PARSE START ====" << std::endl;
#ifdef NDEBUG
  pegtl::dummy_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#else  // NDEBUG
  pegtl::basic_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#endif // NDEBUG
  //std::cout << "==== PARSE END ====" << std::endl << std::endl;

  // Filter out repeated statistics
  unique(std::get<control::STATISTICS>(stack));

  // Store off parsed bundles
  m_control->set(stack);
  m_control->set(boolstack);
}

void
QuinoaParser::unique(std::vector<control::Product>& statistics)
//******************************************************************************
//  Make requested statistics unique
//! \param[in,out]  statistics  Vector of statistics
//! \author  J. Bakosi
//******************************************************************************
{
  std::sort(statistics.begin(), statistics.end());
  auto it = std::unique(statistics.begin(), statistics.end());
  statistics.resize( std::distance(statistics.begin(), it) );
}

void
QuinoaParser::echoGeometry() const
//******************************************************************************
//  Echo parsed data specific to geometry definition
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Geometry: "
            << grammar::Geometry.name(m_control->get<control::GEOMETRY>())
            << std::endl;

  m_control->echoVec<control::BOXES>("Boxes");
}

void
QuinoaParser::echoPhysics() const
//******************************************************************************
//  Echo parsed data specific to physics
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Physics: "
            << grammar::Physics.name(m_control->get<control::PHYSICS>())
            << std::endl;

  m_control->echo<control::NSTEP>("Number of time steps");
  m_control->echo<control::TERM>("Terminate time");
  m_control->echo<control::DT>("Time step size");
  m_control->echo<control::NPAR>("Number of particles");
  m_control->echo<control::TTYI>("TTY output interval");
  m_control->echo<control::DMPI>("Dump output interval");
  m_control->echo<control::STAI>("Statistics output interval");
  m_control->echo<control::PDFI>("PDF output interval");
  m_control->echo<control::GLBI>("Glob output interval");
  m_control->echoVecVecNames<control::STATISTICS>("Requested statistics",true);
  m_control->echoVecVecNames<control::STATISTICS>("Estimated statistics");
  m_control->echo<control::INPUT>("Input filename");
  m_control->echo<control::OUTPUT>("Output filename");
  m_control->echo<control::PDFNAME>("PDF filename");
  m_control->echo<control::GLOBNAME>("Glob filename");
  m_control->echo<control::STATNAME>("Statistics filename");
}

void
QuinoaParser::echoMass() const
//******************************************************************************
//  Echo parsed data specific to mass model
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Mass model: "
            << grammar::Mass.name(m_control->get<control::MASS>())
            << std::endl;

  m_control->echo<control::NDENSITY>("Number of density components");
  m_control->echo<control::AT>("At");
}

void
QuinoaParser::echoHydro() const
//******************************************************************************
//  Echo parsed data specific to hydrodynamics model
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Hydrodynamics model: "
            << grammar::Hydro.name(m_control->get<control::HYDRO>())
            << std::endl;

  m_control->echo<control::NVELOCITY>("Number of velocity components");
  m_control->echo<control::C0>("C0");
}

void
QuinoaParser::echoMix() const
//******************************************************************************
//  Echo parsed data specific to mix model
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Material mix model: "
            << grammar::Mix.name(m_control->get<control::MIX>())
            << std::endl;

  m_control->echo<control::NSCALAR>("Number of scalar components");
  m_control->echoVec<control::B>("Parameter vector b");
  m_control->echoVec<control::S>("Parameter vector S");
  m_control->echoVec<control::KAPPA>("Parameter vector kappa");
  m_control->echoVec<control::C>("Parameter vector c");
}

void
QuinoaParser::echoFrequency() const
//******************************************************************************
//  Echo parsed data specific to turbulence frequency model
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * Turbulence frequency model: "
            << grammar::Frequency.name(m_control->get<control::FREQUENCY>())
            << std::endl;

  m_control->echo<control::NFREQUENCY>("Number of turbulence frequency "
                                       " components");
  m_control->echo<control::FREQ_GAMMA_C1>("C1");
  m_control->echo<control::FREQ_GAMMA_C2>("C2");
  m_control->echo<control::FREQ_GAMMA_C3>("C3");
  m_control->echo<control::FREQ_GAMMA_C4>("C4");
}

void
QuinoaParser::echo() const
//******************************************************************************
//  Echo parsed information from quinoa control
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "Parsed from " << m_filename << ":\n" << std::setfill('-')
            << std::setw(13+m_filename.length()) << "-" << std::endl;

  if (m_control->set<control::TITLE>())
    std::cout << " * Title: " << m_control->get<control::TITLE>() << std::endl;

  if (m_control->set<control::GEOMETRY>()) echoGeometry();
  if (m_control->set<control::PHYSICS>()) echoPhysics();
  if (m_control->set<control::MASS>()) echoMass();
  if (m_control->set<control::HYDRO>()) echoHydro();
  if (m_control->set<control::MIX>()) echoMix();
  if (m_control->set<control::FREQUENCY>()) echoFrequency();

  std::cout << std::endl;
}
