//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Tue 08 Apr 2014 09:30:06 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <Driver.h>
#include <Print.h>
#include <MeshConv/CmdLine/CmdLine.h>

namespace meshconv {

//! MeshConvDriver : Driver
class MeshConvDriver : public tk::Driver {

  public:
    //! Constructor
    explicit MeshConvDriver(int argc, char** argv, const tk::Print& print);

    //! Destructor
    ~MeshConvDriver() noexcept override = default;

    //! Execute
    void execute() const override;

  private:
    //! Don't permit copy constructor
    MeshConvDriver(const MeshConvDriver&) = delete;
    //! Don't permit assigment constructor
    MeshConvDriver& operator=(const MeshConvDriver&) = delete;
    //! Don't permit move constructor
    MeshConvDriver(MeshConvDriver&&) = delete;
    //! Don't permit move assignment
    MeshConvDriver& operator=(MeshConvDriver&&) = delete;

    std::unique_ptr< ctr::CmdLine > m_cmdline;
};

} // meshconv::

#endif // MeshConvDriver_h
