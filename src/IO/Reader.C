//******************************************************************************
/*!
  \file      src/IO/Reader.C
  \author    J. Bakosi
  \date      Mon Oct  7 08:25:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Reader class definition
  \details   Reader class definition
*/
//******************************************************************************

#include <iostream>

#include <Reader.h>
#include <Exception.h>

using namespace tk;

Reader::Reader(const std::string filename) :
  m_filename(filename)
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  // Check if file exists, throw exception if it does not
  m_inFile.open(filename, std::ifstream::in);
  ErrChk(m_inFile.good(), ExceptType::FATAL,
         "Failed to open file: " + filename);

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  m_inFile.get();
  ErrChk(m_inFile.good(), ExceptType::FATAL,
         "Failed to read from file: " + filename);

  // Close it
  m_inFile.close();
  ErrChk(!m_inFile.fail(), ExceptType::FATAL,
         "Failed to close file: " + filename);

  // Re-open
  m_inFile.open(filename, std::ifstream::in);
  ErrChk(m_inFile.good(), ExceptType::FATAL,
         "Failed to open file: " + filename);
}

Reader::~Reader() noexcept
//******************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_inFile.close();
    ErrChk(!m_inFile.fail(), ExceptType::WARNING,
           "Failed to close file: " + m_filename);

  } // emit only a warning on error
    catch (Exception& e) {
      e.echo("WARNING");
    }
    catch (std::exception& e) {
      std::cout << ">>> std::exception in MeshReader destructor: " << e.what()
                << std::endl;
    }
    catch (...) {
      std::cout << ">>> UNKNOWN EXCEPTION in MeshReader destructor"
                << std::endl;
    }
}
