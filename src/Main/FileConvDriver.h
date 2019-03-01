// *****************************************************************************
/*!
  \file      src/Main/FileConvDriver.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     File converter driver
  \details   File converter driver.
*/
// *****************************************************************************
#ifndef FileConvDriver_h
#define FileConvDriver_h

#include <iosfwd>

#include "FileConv/CmdLine/CmdLine.h"

namespace tk { class Print; }

//! File converter declarations and definitions
namespace fileconv {

//! File Converter driver used polymorphically with tk::Driver
class FileConvDriver {

  public:
    //! Constructor
    explicit FileConvDriver( const tk::Print&,
                             const ctr::CmdLine& cmdline );

    //! Execute
    void execute() const;

  private:

    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // fileconv::

#endif // FileConvDriver_h
