//******************************************************************************
/*!
  \file      src/Main/MeshConv.C
  \author    J. Bakosi
  \date      Tue 26 Aug 2014 12:19:14 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gmsh to Exodus II mesh file converter
  \details   Gmsh to Exodus II mesh file converter
*/
//******************************************************************************

#include <Config.h>
#include <MeshConvDriver.h>
#include <TPLInfo/ExodusII.h>
#include <TPLInfo/MKL.h>
#include <TPLInfo/Boost.h>
#include <MeshConv/CmdLine/Parser.h>
#include <meshconv.decl.h>
#include <Init.h>

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
#ifdef HAS_MKL
  echoMKL( print, "Intel Math Kernel Library" );
#else
  print.item( "Intel Math Kernel Library", "n/a" );
#endif
  echoBoost( print, "Boost C++ Libraries" );
  echoExodusII( print, "ExodusII library" );
  print.endpart();
}

} // meshconv::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg )
    try :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
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
    } catch (...) { processException(); }

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

    //! Process an exception
    void processException() {
      try {
        throw;      // rethrow exception to deal with it here
      }
        // Catch Quina::Exceptions
        catch ( tk::Exception& qe ) {
          qe.handleException();
        }
        // Catch std::exception and transform it into Quinoa::Exception without
        // file:line:func information
        catch ( std::exception& se ) {
          tk::Exception qe( se.what() );
          qe.handleException();
        }
        // Catch uncaught exception
        catch (...) {
          tk::Exception qe( "Non-standard exception" );
          qe.handleException();
        }

      // Tell the runtime system to exit
      finalize();
    }

  private:
    meshconv::ctr::CmdLine m_cmdline;           //!< Command line
    meshconv::CmdLineParser m_cmdParser;        //!< Command line parser
    tk::Print m_print;                          //!< Pretty printer
    meshconv::MeshConvDriver m_driver;          //!< Driver
    std::vector< tk::Timer > m_timer;           //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <meshconv.def.h>
