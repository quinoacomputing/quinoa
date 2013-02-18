//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 02:26:35 PM MST
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
#include <IOException.h>
#include <ParserException.h>

using namespace Quinoa;
using namespace control;

Parser::Parser(const string& filename, Control* const control) :
  m_filename(filename), m_control(control)
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \param[in]  control       Control category where parsed data are stored
//! \author  J. Bakosi
//******************************************************************************
{
  // Check if control file exists, throw exception if not
  m_q.open(m_filename, ifstream::in);
  Assert(m_q.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
  m_q.close();
  Assert(!m_q.fail(), IOException,FATAL,IO_FAILED_CLOSE,m_filename);
}

void
Parser::parse()
//******************************************************************************
//  Parse
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize new bundle for sparsed data with defaults
  Bundle stack(Defaults);

  //cout << "==== PARSE START ====" << endl;
#ifdef NDEBUG
  pegtl::dummy_parse_file< grammar::read_file >( m_filename, stack );
#else  // NDEBUG
  pegtl::basic_parse_file< grammar::read_file >( m_filename, stack );
#endif // NDEBUG
  //cout << "==== PARSE END ====" << endl << endl;

  // Store off stuff parsed
  m_control->set(stack);
}

void
Parser::echo()
//******************************************************************************
//  Echo information on stuff parsed
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Parsed from " << m_filename << ":\n" << setfill('-')
       << setw(13+m_filename.length()) << "-" << endl;

  if (m_control->set<TITLE>()) {
    cout << " * Title: " << m_control->get<TITLE>() << endl;
  }

  if (m_control->set<PHYSICS>()) {
    cout << " * Physics: " << m_control->physicsName() << endl;
  }

  if (m_control->set<HYDRO>()) {
    cout << " * Hydrodynamics: " << m_control->hydroName() << endl;
  }

  if (m_control->set<MIX>()) {

    cout << " * Material mixing: " << m_control->mixName() << endl;

    if (m_control->set<NSTEP>()) {
      cout << "   - Number of time steps: " << m_control->get<NSTEP>() << endl;
    }

    if (m_control->set<TERM>()) {
      cout << "   - Terminate time: " << m_control->get<TERM>() << endl;
    }

    if (m_control->set<DT>()) {
      cout << "   - Time step size: " << m_control->get<DT>() << endl;
    }

    if (m_control->set<NSCALAR>()) {
      cout << "   - Number of mixing scalars: " << m_control->get<NSCALAR>()
           << endl;
    }

    if (m_control->set<NPAR>()) {
      cout << "   - Number of particles: " << m_control->get<NPAR>() << endl;
    }

    if (m_control->set<ECHO>()) {
      cout << "   - Screen-output every " << m_control->get<ECHO>() << " step"
           << endl;
    }

  }

  cout << endl;
}
