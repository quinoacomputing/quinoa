//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Tue May  7 17:17:06 2013
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

Parser::Parser(const string& filename, Control* const control) :
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
  // Check if control file exists, throw exception if not
  m_q.open(m_filename, ifstream::in);
  ErrChk(m_q.good(), FATAL, "Failed to open file: " + m_filename);

  m_q.close();
  ErrChk(!m_q.fail(), FATAL, "Failed to close file: " + m_filename);
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
    boolstack(tuple_size<decltype(control::DEFAULTS)>::value, false);

  //cout << "==== PARSE START ====" << endl;
#ifdef NDEBUG
  pegtl::dummy_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#else  // NDEBUG
  pegtl::basic_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#endif // NDEBUG
  //cout << "==== PARSE END ====" << endl << endl;

  // Filter out repeated statistics
  unique(get<control::STATISTICS>(stack));

  // Store off parsed bundles
  m_control->set(stack);
  m_control->set(boolstack);
}

void
Parser::unique(vector<control::Product>& statistics)
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
Parser::echoMix() const
//******************************************************************************
//  Echo parsed data from material mix block
//! \author  J. Bakosi
//******************************************************************************
{
  cout << " * Material mixing: " << m_control->mixName() << endl;

  m_control->echo<control::NSTEP>("Number of time steps");
  m_control->echo<control::TERM>("Terminate time");
  m_control->echo<control::DT>("Time step size");
  m_control->echo<control::NSCALAR>("Number of mixing scalars");
  m_control->echo<control::NPAR>("Number of particles");
  m_control->echo<control::TTYI>("TTY output interval");
  m_control->echo<control::DUMP>("Dump output interval");
  m_control->echo<control::PLTI>("Plot output interval");
  m_control->echo<control::PDFI>("PDF output interval");
  m_control->echo<control::GLOB>("Glob output interval");
  m_control->echo<control::C0>("C0");
  m_control->echoVec<control::B>("Parameter vector b");
  m_control->echoVec<control::S>("Parameter vector S");
  m_control->echoVec<control::KAPPA>("Parameter vector kappa");
  m_control->echoVec<control::C>("Parameter vector c");
  m_control->echoVecVecNames<control::STATISTICS>("Requested statistics",true);
  m_control->echoVecVecNames<control::STATISTICS>("Estimated statistics");
  m_control->echo<control::PDFNAME>("PDF base filename");
  m_control->echo<control::GLOBNAME>("Glob filename");
  m_control->echo<control::PLOTNAME>("Plot base filename");
}

void
Parser::echoHydro() const
//******************************************************************************
//  Echo parsed data from hydrodynamics block
//! \author  J. Bakosi
//******************************************************************************
{
  cout << " * Hydrodynamics: " << m_control->hydroName() << endl;

  m_control->echo<control::NSTEP>("Number of time steps");
  m_control->echo<control::TERM>("Terminate time");
  m_control->echo<control::DT>("Time step size");
  m_control->echo<control::NPAR>("Number of particles");
  m_control->echo<control::TTYI>("TTY output interval");
  m_control->echo<control::DUMP>("Dump output interval");
  m_control->echo<control::PLTI>("Plot output interval");
  m_control->echo<control::PDFI>("PDF output interval");
  m_control->echo<control::GLOB>("Glob output interval");
  m_control->echo<control::C0>("C0");
  m_control->echoVecVecNames<control::STATISTICS>("Requested statistics",true);
  m_control->echoVecVecNames<control::STATISTICS>("Estimated statistics");
  m_control->echo<control::PDFNAME>("PDF base filename");
  m_control->echo<control::GLOBNAME>("Glob filename");
  m_control->echo<control::PLOTNAME>("Plot base filename");
}

void
Parser::echo() const
//******************************************************************************
//  Echo information on stuff parsed
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Parsed from " << m_filename << ":\n" << setfill('-')
       << setw(13+m_filename.length()) << "-" << endl;

  if (m_control->set<control::TITLE>())
    cout << " * Title: " << m_control->get<control::TITLE>() << endl;

  if (m_control->set<control::PHYSICS>())
    cout << " * Physics: " << m_control->physicsName() << endl;

  if (m_control->set<control::HYDRO>()) echoHydro();

  if (m_control->set<control::MIX>()) echoMix();

  cout << endl;
}
