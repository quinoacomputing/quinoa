//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Sun 08 Sep 2013 09:10:53 PM MDT
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
QuinoaParser::echo() const
//******************************************************************************
//  Echo parsed information from quinoa control
//! \author  J. Bakosi
//******************************************************************************
{
  using namespace control;
  std::cout << "Parsed from " << m_filename << ":\n" << std::setfill('-')
            << std::setw(13+m_filename.length()) << "-" << std::endl;

  std::cout << " * Title: " << m_control.get<title>() << std::endl;
  std::cout << " * Geometry: "
            << grammar::Geometry.name(m_control.get<selected,geometry>())
            << std::endl;
  std::cout << " * Physics: "
            << grammar::Physics.name(m_control.get<selected,physics>())
            << std::endl;
  std::cout << " * Position: "
            << grammar::Mass.name(m_control.get<selected,mass>())
            << std::endl;
  std::cout << " * Hydrodynamics: "
            << grammar::Hydro.name(m_control.get<selected,hydro>())
            << std::endl;
  std::cout << " * Energy: "
            << grammar::Energy.name(m_control.get<selected,energy>())
            << std::endl;
  std::cout << " * Material mix: "
            << grammar::Mix.name(m_control.get<selected,mix>())
            << std::endl;
  std::cout << " * Frequency: "
            << grammar::Frequency.name(m_control.get<selected,frequency>())
            << std::endl;
  std::cout << " * Material mix rate: "
            << grammar::MixRate.name(m_control.get<selected,mixrate>())
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

  std::cout << std::endl;
}
