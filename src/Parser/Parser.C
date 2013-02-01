//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Fri 01 Feb 2013 06:06:39 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Control.h>
#include <Grammar.h>
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
  cout << "==== PARSE START ====" << endl;
  grammar::stack_type stack;
#ifdef NDEBUG
  pegtl::dummy_parse_file< grammar::read_file >( m_filename, stack );
#else  // NDEBUG
  pegtl::basic_parse_file< grammar::read_file >( m_filename, stack );
#endif // NDEBUG
  cout << "==== PARSE END ====" << endl << endl;

  // Store off stuff parsed
  m_control->setTitle(stack.title);
  m_control->setPhysics(stack.physics);
  m_control->setHydro(stack.hydro);
  m_control->setMix(stack.mix);

  cout << m_control->title() << endl;
  cout << static_cast<int>(m_control->physics()) << endl;
  cout << static_cast<int>(m_control->hydro()) << endl;
  cout << static_cast<int>(m_control->mix()) << endl;
}

void
Parser::echo()
//******************************************************************************
//  Echo information on stuff parsed
//! \author  J. Bakosi
//******************************************************************************
{
}
