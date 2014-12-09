//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:31:07 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's command line
  \details   Quinoa's command line
*/
//******************************************************************************
#ifndef QuinoaCmdLine_h
#define QuinoaCmdLine_h

#include <string>

#include <Control.h>
#include <Quinoa/Types.h>

namespace quinoa {
namespace ctr {

//! CmdLine : Control< specialized to Quinoa >, see Types.h,
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
      set< tag::io, tag::glob >( "glob.txt" );
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
} // quinoa::

#endif // QuinoaCmdLine_h
