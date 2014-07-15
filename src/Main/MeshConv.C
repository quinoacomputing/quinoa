//******************************************************************************
/*!
  \file      src/Main/MeshConv.C
  \author    J. Bakosi
  \date      Tue 15 Jul 2014 08:59:17 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Gmsh to Exodus II mesh file converter
  \details   Gmsh to Exodus II mesh file converter
*/
//******************************************************************************

#include <Config.h>
#include <Init.h>
#include <MeshConvDriver.h>
#include <TPLInfo/ExodusII.h>
#include <MeshConv/CmdLine/Parser.h>
#include <meshconv.decl.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

namespace meshconv {

void echoTPL( const tk::Print& print )
//******************************************************************************
//  Echo TPL version informaion for libs specific to MeshConv
//! \author  J. Bakosi
//******************************************************************************
{
  print << '\n';
  echoExodusII( print, "ExodusII library" );
}

} // meshconv::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg ) :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : tk::null ),
      // Create MeshConv driver
      m_driver( tk::Main< meshconv::MeshConvDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::MESHCONV,
                          MESHCONV_EXECUTABLE,
                          m_print,
                          meshconv::echoTPL ) ),
      m_timer(1)        // Start new timer measuring the total runtime
    {
      delete msg;
      mainProxy = thisProxy;
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
      // Start new timer measuring the migration of global-scope data
      m_timer.emplace_back();
    }

    void execute() {
      m_timestamp.emplace( "Migration of global-scope data", m_timer[1].hms() );
      m_driver.execute();       // does not fire up async chares
      finalize();
    }

    void finalize() {
      m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
      m_print.time( "Timers (h:m:s)", m_timestamp );
      m_print.endpart();
      CkExit();
    }

  private:
    meshconv::ctr::CmdLine m_cmdline;           //!< Command line
    meshconv::CmdLineParser m_cmdParser;        //!< Command line parser
    tk::Print m_print;                          //!< Pretty printer
    meshconv::MeshConvDriver m_driver;          //!< Driver
    std::vector< tk::Timer > m_timer;           //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <meshconv.def.h>
