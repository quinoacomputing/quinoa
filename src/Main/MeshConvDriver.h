//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Sun 31 May 2015 06:24:58 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <iosfwd>

#include "MeshConv/CmdLine/CmdLine.h"

namespace tk { class Print; }

//! Mesh converter declarations and definitions
namespace meshconv {

//! Mesh converter driver used polymorphically with tk::Driver
class MeshConvDriver {

  public:
    //! Constructor
    explicit MeshConvDriver( const tk::Print& print,
                             const ctr::CmdLine& cmdline );

    //! Execute
    void execute() const;

  private:
    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // meshconv::

#endif // MeshConvDriver_h
