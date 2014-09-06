//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Fri 05 Sep 2014 07:55:48 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Text statistics writer
  \details   Text statistics writer
*/
//******************************************************************************
#ifndef TxtStatWriter_h
#define TxtStatWriter_h

#include <string>
#include <vector>

#include <Types.h>
#include <Writer.h>

namespace quinoa {

//! TxtStatWriter : Writer
class TxtStatWriter : public tk::Writer {

  public:
    //! Constructor
    explicit TxtStatWriter( const std::string& filename ) :
      Writer( filename ) {}

    //! Write out statistics file header
    void header( const std::vector< bool >& plotOrd,
                 const std::vector< std::string >& nameOrd,
                 const std::vector< std::string >& nameCen ) const;

    //! Write statistics file
    void writeStat( int it,
                    tk::real t,
                    const std::vector< tk::real >& ordinary,
                    const std::vector< tk::real >& central,
                    const std::vector< bool >& plotOrd );
};

} // quinoa::

#endif // TxtStatWriter_h
