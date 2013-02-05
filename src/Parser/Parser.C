//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Mon 04 Feb 2013 09:43:05 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Grammar.h>
#include <Control.h>
#include <Parser.h>
#include <IOException.h>
#include <ParserException.h>

using namespace Quinoa;

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
  //cout << "==== PARSE START ====" << endl;
  control::Bundle stack;
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

  cout << " * Title: " << m_control->get<control::TITLE>() << endl;

  if (m_control->get<control::PHYSICS>() > control::NO_PHYSICS) {
    cout << " * Physics: " << m_control->physicsName() << endl;
  }

  if (m_control->get<control::HYDRO>() > control::NO_HYDRO) {
    cout << " * Hydrodynamics: " << m_control->hydroName() << endl;
  }

  if (m_control->get<control::MIX>() > control::NO_MIX) {
    cout << " * Material mixing: " << m_control->mixName() << endl;
    cout << "   - Number of time steps: " << m_control->get<control::NSTEP>()
         << endl;
    cout << "   - Terminate at time: " << m_control->get<control::TERM>()
         << endl;
    cout << "   - Size of time step: " << m_control->get<control::DT>()
         << endl;
    cout << "   - Number of mixing scalars: "
         << m_control->get<control::NSCALAR>() << endl;
    cout << "   - Number of particles: " << m_control->get<control::NPAR>()
         << endl;
    cout << "   - Screen-output at every " << m_control->get<control::ECHO>()
         << " step" << endl;
  }

  cout << endl;
}
