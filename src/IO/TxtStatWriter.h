//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 03:56:06 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Text statistics writer
  \details   Text statistics writer
*/
//******************************************************************************
#ifndef TxtStatWriter_h
#define TxtStatWriter_h

#include <string>
#include <vector>
#include <fstream>

#include <Types.h>
#include <Writer.h>

namespace tk {

//! TxtStatWriter : Writer
class TxtStatWriter : public tk::Writer {

  public:
    //! Constructor
    explicit TxtStatWriter( const std::string& filename,
                            std::ios_base::openmode mode = std::ios_base::out )
      : Writer( filename, mode ) {}

    //! Write out statistics file header
    void header( const std::vector< std::string >& nameOrd,
                 const std::vector< std::string >& nameCen ) const;

    //! Write statistics file
    std::size_t stat( int it, tk::real t,
                      const std::vector< tk::real >& ordinary,
                      const std::vector< tk::real >& central );
};

} // tk::

#endif // TxtStatWriter_h
