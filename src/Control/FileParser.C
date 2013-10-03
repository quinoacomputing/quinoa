//******************************************************************************
/*!
  \file      src/Control/FileParser.C
  \author    J. Bakosi
  \date      Thu Oct  3 12:01:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************

#include <fstream>

#include <FileParser.h>
#include <Exception.h>

using namespace quinoa;

FileParser::FileParser(Base& base, const std::string& filename) :
  Parser(base),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[inout]  base          Essentials
//! \param[in]     filename      File to parse
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  std::ifstream q;

  // Check if file exists, throw exception if it does not
  q.open(m_filename, std::ifstream::in);
  ErrChk(q.good(), ExceptType::FATAL, "Failed to open file: " + m_filename);

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the  read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  q.get();
  ErrChk(q.good(), ExceptType::FATAL, "Failed to read from file: " + m_filename);

  // Close it
  q.close();
  ErrChk(!q.fail(), ExceptType::FATAL, "Failed to close file: " + m_filename);
}
