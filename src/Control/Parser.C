//******************************************************************************
/*!
  \file      src/Control/Parser.C
  \author    J. Bakosi
  \date      Thu 26 Sep 2013 08:40:21 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <fstream>

#include <Parser.h>
#include <Exception.h>

using namespace quinoa;

Parser::Parser(const std::string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  std::ifstream q;

  // Check if control file exists, throw exception if it does not
  q.open(m_filename, std::ifstream::in);
  ErrChk(q.good(), ExceptType::FATAL, "Failed to open file: " + m_filename);

  q.close();
  ErrChk(!q.fail(), ExceptType::FATAL, "Failed to close file: " + m_filename);
}
