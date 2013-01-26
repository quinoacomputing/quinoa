//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 10:03:28 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Grammar.def.h>
#include <Parser.h>
#include <Control.h>
#include <IOException.h>

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
//  m_q.open(m_filename, ifstream::in);
//  Assert(m_q.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
}

Parser::~Parser()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
//  m_q.close();

  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
//  if (m_q.fail())
//    cout << "WARNING: Failed to close file: " << m_filename << endl;
}

void
Parser::Parse()
//******************************************************************************
//  Parse
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "==== PARSE START ====" << endl;
  pegtl::basic_parse_file< grammar::read_file >( m_filename );
  cout << "==== PARSE END ====" << endl << endl;

  m_control->setTitle( grammar::title );

  cout << "Title: " << m_control->title() << endl << endl;
}
