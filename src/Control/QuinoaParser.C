//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Sun 08 Sep 2013 03:31:42 PM MDT
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
  using namespace control;

  //std::cout << "==== PARSE START ====" << std::endl;
#ifdef NDEBUG
  //pegtl::dummy_parse_file<grammar::read_file>(m_filename, m_control);
#else  // NDEBUG
  //pegtl::basic_parse_file<grammar::read_file>(m_filename, m_control);
#endif // NDEBUG
  //std::cout << "==== PARSE END ====" << std::endl << std::endl;

  // Filter out repeated statistics
  unique(const_cast<std::vector<Product>&>(m_control.get<stats>()));
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
  statistics.resize(std::distance(statistics.begin(), it));
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
            << grammar::Geometry.name(m_control.get<selected,geometry>())
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
            << grammar::Physics.name(m_control.get<selected,physics>())
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
  m_control.echoVecVecNames<stats>("Requested statistics",true);
  m_control.echoVecVecNames<stats>("Estimated statistics");
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
            << grammar::Mass.name(m_control.get<selected,mass>())
            << std::endl;

  m_control.echo<component,ndensity>("Number of density components");
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
            << grammar::Hydro.name(m_control.get<selected,hydro>())
            << std::endl;

  m_control.echo<component,nvelocity>("Number of velocity components");
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
            << grammar::Mix.name(m_control.get<selected,mix>())
            << std::endl;

  m_control.echo<component,nscalar>("Number of scalar components");
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
            << grammar::Frequency.name(m_control.get<selected,frequency>())
            << std::endl;

  m_control.echo<component,nfrequency>("Number of turbulence frequency "
                                       " components");
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

  m_control.echo<title>(" * Title: ");

  echoGeometry();
  echoPhysics();
  echoMass();
  echoHydro();
  echoMix();
  echoFrequency();

  std::cout << std::endl;
}
