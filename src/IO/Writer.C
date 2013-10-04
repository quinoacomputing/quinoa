//******************************************************************************
/*!
  \file      src/IO/Writer.C
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 08:30:17 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Writer base class definition
  \details   Writer base class definition
*/
//******************************************************************************

#include <iostream>

#include <Writer.h>
#include <Exception.h>

using namespace quinoa;

Writer::Writer(const std::string& filename) :
  m_filename(filename),
  m_outFile()
//******************************************************************************
//  Constructor: Acquire file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile.open(m_filename, std::ofstream::out);
  ErrChk(m_outFile.good(), ExceptType::FATAL,
         "Failed to open file: " + m_filename);
}

Writer::~Writer() noexcept
//******************************************************************************
//  Destructor: Release file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_outFile.close();

    if (m_outFile.fail())
      std::cout << "WARNING: Failed to close file: " << m_filename << std::endl;

  } // emit only a warning on error
    catch (Exception& e) {
      e.echo("WARNING");
    }
    catch (std::exception& e) {
      std::cout << ">>> std::exception in Writer destructor: " << e.what()
                << std::endl;
    }
    catch (...) {
      std::cout << "UNKNOWN EXCEPTION in Writer destructor" << std::endl;
    }
}
