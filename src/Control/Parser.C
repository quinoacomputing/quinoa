//******************************************************************************
/*!
  \file      src/Control/Parser.C
  \author    J. Bakosi
  \date      Wed Aug 28 15:20:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <fstream>

#include <Parser.h>
#include <Exception.h>

using namespace Quinoa;

Parser::Parser(const std::string& filename, Control* const control) :
  m_filename(filename),
  m_control(control)
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \param[in]  control       Control category where parsed data are stored
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
