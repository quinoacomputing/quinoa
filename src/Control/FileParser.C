//******************************************************************************
/*!
  \file      src/Control/FileParser.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 06:34:35 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************

#include <fstream>

#include <FileParser.h>
#include <Exception.h>

using tk::FileParser;

FileParser::FileParser(const std::string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]     filename      File to parse
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream is in a valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  //! Make sure there is a filename
  Assert(!filename.empty(), ExceptType::FATAL, "No filename specified");

  std::ifstream q;

  // Check if file exists, throw exception if it does not
  q.open(filename, std::ifstream::in);
  ErrChk(q.good(), ExceptType::FATAL, "Failed to open file: " + filename);

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  q.get();
  ErrChk(q.good(), ExceptType::FATAL, "Failed to read from file: " + filename);

  // Close it
  q.close();
  ErrChk(!q.fail(), ExceptType::FATAL, "Failed to close file: " + filename);
}
