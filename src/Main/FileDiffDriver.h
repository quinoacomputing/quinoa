// *****************************************************************************
/*!
  \file      src/Main/FileDiffDriver.h
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     File converter driver
  \details   File converter driver.
*/
// *****************************************************************************
#ifndef FileDiffDriver_h
#define FileDiffDriver_h

#include <iosfwd>

#include "FileDiff/CmdLine/CmdLine.h"

namespace tk { class Print; }

//! File converter declarations and definitions
namespace filediff {

//! File Difference driver used polymorphically with tk::Driver
class FileDiffDriver {

  public:
    //! Constructor
    explicit FileDiffDriver( const tk::Print& print,
                             const ctr::CmdLine& cmdline );

    //! Execute
    void execute() const;

  private:

    const tk::Print& m_print;           //!< Pretty printer
    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // filediff::

#endif // FileDiffDriver_h
