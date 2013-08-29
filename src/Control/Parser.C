//******************************************************************************
/*!
  \file      src/Control/Parser.C
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 07:48:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <fstream>

#include <Parser.h>
#include <Exception.h>

using namespace Quinoa;

Parser::Parser(const std::string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  std::ifstream m_q;

  // Check if control file exists, throw exception if it does not
  m_q.open(m_filename, std::ifstream::in);
  ErrChk(m_q.good(), ExceptType::FATAL, "Failed to open file: " + m_filename);

  m_q.close();
  ErrChk(!m_q.fail(), ExceptType::FATAL, "Failed to close file: " + m_filename);
}
