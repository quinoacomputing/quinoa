//******************************************************************************
/*!
  \file      src/IO/Reader.C
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 09:50:51 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief      reader class definition
  \details    reader class definition
*/
//******************************************************************************

#include <iostream>

#include <Reader.h>

using namespace Quinoa;

Reader::Reader(const std::string filename) :
  m_filename(filename),
  m_inFile()
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inFile.open(m_filename, std::ifstream::in);
  ErrChk(m_inFile.good(), ExceptType::FATAL,
         "Failed to open file: " + m_filename);
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
