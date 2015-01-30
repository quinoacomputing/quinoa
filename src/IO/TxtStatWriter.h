//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 10:14:24 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Text statistics writer declaration
  \details   This file declares the ASCII statistics writer class that
     facilitates outputing statistics to text files.
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

//! \brief TxtStatWriter : tk::Writer
//! \details ASCII statistics writer class that facilitates outputing
//!   statistics to text files.
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
