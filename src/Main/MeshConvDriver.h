// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************
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
    const tk::Print& m_print;           //!< Pretty printer
    const bool m_reorder;               //!< Whether to also reorder mesh nodes
    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // meshconv::

#endif // MeshConvDriver_h
