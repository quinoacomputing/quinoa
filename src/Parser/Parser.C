//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Tue 30 Jul 2013 08:11:03 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Grammar.h>
#include <Control.h>
#include <Defaults.h>
#include <Parser.h>
#include <Exception.h>

using namespace Quinoa;

Parser::Parser(const std::string& filename, Control* const control) :
  m_filename(filename),
  m_control(control),
  m_q()
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \param[in]  control       Control category where parsed data are stored
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  // Check if control file exists, throw exception if it does not
  m_q.open(m_filename, std::ifstream::in);
  ErrChk(m_q.good(), ExceptType::FATAL, "Failed to open file: " + m_filename);

  m_q.close();
  ErrChk(!m_q.fail(), ExceptType::FATAL, "Failed to close file: " + m_filename);
}

void
Parser::parse()
//******************************************************************************
//  Parse
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize new bundle for parsed data with defaults
  control::Bundle stack(control::DEFAULTS);
  // Initialize new bool bundle for indicating what data is set in bundle
  control::BoolBundle
    boolstack(std::tuple_size<decltype(control::DEFAULTS)>::value, false);

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
Parser::unique(std::vector<control::Product>& statistics)
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
Parser::echoGeometry() const
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
Parser::echoPhysics() const
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
Parser::echoMass() const
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
Parser::echoHydro() const
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
Parser::echoMix() const
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
Parser::echoFrequency() const
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
Parser::echoRNGTest() const
//******************************************************************************
//  Echo parsed data specific to RNG test suite
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << " * RNG test suite: "
            << grammar::RNGTest.name(m_control->get<control::RNGTEST>())
            << std::endl;
}

void
Parser::echo() const
//******************************************************************************
//  Echo information on stuff parsed
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
  if (m_control->set<control::RNGTEST>()) echoRNGTest();

  std::cout << std::endl;
}
