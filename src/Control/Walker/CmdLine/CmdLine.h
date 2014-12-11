//******************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:28:47 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's command line
  \details   Walker's command line
*/
//******************************************************************************
#ifndef WalkerCmdLine_h
#define WalkerCmdLine_h

#include <string>

#include <Control.h>
#include <Walker/Types.h>

namespace walker {
//! Walker control facilitating user input to internal data transfer
namespace ctr {

//! CmdLine : Control< specialized to Walker >, see Types.h,
class CmdLine : public tk::Control<
                  // tag               type
                  tag::io,             ios,
                  tag::virtualization, tk::real,
                  tag::verbose,        bool,
                  tag::error,          std::vector< std::string > > {

  public:
    //! Constructor: set all defaults. Anything not set here initialized by the
    //! compiler using the default constructor for the corresponding type.
    CmdLine() {
      set< tag::io, tag::output >( "out" );
      set< tag::io, tag::pdf >( "pdf" );
      set< tag::io, tag::stat >( "stat.txt" );
      set< tag::virtualization >( 0.0 );
      set< tag::verbose >( false ); // Quiet output by default
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::io,             ios,
                   tag::virtualization, tk::real,
                   tag::verbose,        bool,
                   tag::error,          std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
};

} // ctr::
} // walker::

#endif // WalkerCmdLine_h
