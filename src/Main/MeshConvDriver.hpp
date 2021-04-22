// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <iosfwd>

#include "Print.hpp"
#include "MeshConv/CmdLine/CmdLine.hpp"

//! Mesh converter declarations and definitions
namespace meshconv {

//! Mesh converter driver used polymorphically with tk::Driver
class MeshConvDriver {

  public:
    //! Constructor
    explicit MeshConvDriver( const ctr::CmdLine& cmdline, int );

    //! Execute
    void execute() const;

  private:
    const tk::Print m_print;            //!< Pretty printer
    const bool m_reorder;               //!< Whether to also reorder mesh nodes
    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // meshconv::

#endif // MeshConvDriver_h
