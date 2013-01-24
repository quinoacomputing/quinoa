//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Tue 22 Jan 2013 10:42:35 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <iostream>

#include <Parser.h>
#include <IOException.h>

using namespace Quinoa;

Parser::Parser(const string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  m_q.open(m_filename, ifstream::in);
  Assert(m_q.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);  
}

Parser::~Parser()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  m_q.close();

  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_q.fail())
    cout << "WARNING: Failed to close file: " << m_filename << endl;
}
