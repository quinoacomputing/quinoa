//******************************************************************************
/*!
  \file      src/Control/FileParser.C
  \author    J. Bakosi
  \date      Wed Oct  2 15:42:55 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************

#include <fstream>

#include <FileParser.h>
#include <Exception.h>

using namespace quinoa;

FileParser::FileParser(const std::string& filename, Base& base) :
  Parser(base),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]     filename      File to parse
//! \param[inout]  base          Essentials
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  std::ifstream q;

  // Check if file exists, throw exception if it does not
  q.open(m_filename, std::ifstream::in);
  ErrChk(q.good(), ExceptType::FATAL, "Failed to open file: " + m_filename);

  q.close();
  ErrChk(!q.fail(), ExceptType::FATAL, "Failed to close file: " + m_filename);
}
