//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Wed Sep  4 10:04:13 2013
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
  //std::cout << "==== PARSE START ====" << std::endl;
#ifdef NDEBUG
  //pegtl::dummy_parse_file<grammar::read_file>(m_filename, m_control);
#else  // NDEBUG
  //pegtl::basic_parse_file<grammar::read_file>(m_filename, m_control);
#endif // NDEBUG
  //std::cout << "==== PARSE END ====" << std::endl << std::endl;

  // Filter out repeated statistics
  //unique(m_control.get<control::statistic>().get<control::stats>());
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
  using namespace control;
  std::cout << " * Geometry: "
            << grammar::Geometry.name(m_control.get<selected>().get<geometry>())
            << std::endl;
}

void
QuinoaParser::echoPhysics() const
//******************************************************************************
//  Echo parsed data specific to physics
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << " * Physics: "
            << grammar::Physics.name(m_control.get<selected>().get<physics>())
            << std::endl;

  m_control.echo<incpar,nstep>("Number of time steps");
  m_control.echo<incpar,term>("Terminate time");
  m_control.echo<incpar,dt>("Time step size");
  m_control.echo<component,npar>("Number of particles");
  m_control.echo<interval,tty>("TTY output interval");
  m_control.echo<interval,dump>("Dump output interval");
  m_control.echo<interval,plot>("Statistics output interval");
  m_control.echo<interval,pdf>("PDF output interval");
  m_control.echo<interval,glob>("Glob output interval");
  //m_control.echoVecVecNames<control::STATISTICS>("Requested statistics",true);
  //m_control.echoVecVecNames<control::STATISTICS>("Estimated statistics");
  m_control.echo<io,input>("Input filename");
  m_control.echo<io,output>("Output filename");
  m_control.echo<io,pdf>("PDF filename");
  m_control.echo<io,glob>("Glob filename");
  m_control.echo<io,stats>("Statistics filename");
}

void
QuinoaParser::echoMass() const
//******************************************************************************
//  Echo parsed data specific to mass model
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << " * Mass model: "
            << grammar::Mass.name(m_control.get<selected>().get<mass>())
            << std::endl;

  m_control.echo<component,ndensity>("Number of density components");
  //m_control.echo<component,parameter,atwood>("At");
}

void
QuinoaParser::echoHydro() const
//******************************************************************************
//  Echo parsed data specific to hydrodynamics model
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << " * Hydrodynamics model: "
            << grammar::Hydro.name(m_control.get<selected>().get<hydro>())
            << std::endl;

  m_control.echo<component,nvelocity>("Number of velocity components");
  //m_control->echo<control::C0>("C0");
}

void
QuinoaParser::echoMix() const
//******************************************************************************
//  Echo parsed data specific to mix model
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << " * Material mix model: "
            << grammar::Mix.name(m_control.get<selected>().get<mix>())
            << std::endl;

  m_control.echo<component,nscalar>("Number of scalar components");
  //m_control->echoVec<control::B>("Parameter vector b");
  //m_control->echoVec<control::S>("Parameter vector S");
  //m_control->echoVec<control::KAPPA>("Parameter vector kappa");
  //m_control->echoVec<control::C>("Parameter vector c");
}

void
QuinoaParser::echoFrequency() const
//******************************************************************************
//  Echo parsed data specific to turbulence frequency model
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << " * Turbulence frequency model: "
            << grammar::Frequency.name(m_control.get<selected>().get<frequency>())
            << std::endl;

  m_control.echo<component,nfrequency>("Number of turbulence frequency "
                                       " components");
//   m_control->echo<control::FREQ_GAMMA_C1>("C1");
//   m_control->echo<control::FREQ_GAMMA_C2>("C2");
//   m_control->echo<control::FREQ_GAMMA_C3>("C3");
//   m_control->echo<control::FREQ_GAMMA_C4>("C4");
}

void
QuinoaParser::echo() const
//******************************************************************************
//  Echo parsed information from quinoa control
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << "Parsed from " << m_filename << ":\n" << std::setfill('-')
            << std::setw(13+m_filename.length()) << "-" << std::endl;

  //if (m_control->set<control::TITLE>())
    //std::cout << " * Title: " << m_control.get<title>() << std::endl;
    m_control.echo<title>(" * Title: ");

//   if (m_control->set<control::GEOMETRY>()) echoGeometry();
//   if (m_control->set<control::PHYSICS>()) echoPhysics();
//   if (m_control->set<control::MASS>()) echoMass();
//   if (m_control->set<control::HYDRO>()) echoHydro();
//   if (m_control->set<control::MIX>()) echoMix();
//   if (m_control->set<control::FREQUENCY>()) echoFrequency();

  std::cout << std::endl;
}
