// *****************************************************************************
/*!
  \file      src/IO/H5PartWriter.h
  \author    J. Bakosi
  \date      Thu 28 Jul 2016 08:32:55 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     H5Part particles data writer
  \details   H5Part particles data writer class declaration, facilitating
    writing particle coordinates and associated particle fields into HDF5-based
    H5Part data files in parallel, using MPI-IO.
*/
// *****************************************************************************
#ifndef H5PartWriter_h
#define H5PartWriter_h

#include <NoWarning/H5Part.h>

#include "Types.h"

namespace tk {

//! H5Part particles data data writer
//! \details Particles data writer class facilitating writing particle
//!   coordinates and associated particle fields into HDF5-based H5Part data
//!   files in parallel, using MPI-IO.
//! \see http://vis.lbl.gov/Research/H5Part/
class H5PartWriter {

  public:
    //! Constructor: create/open H5Part file
    explicit H5PartWriter( const std::string& filename );

    //! Destructor
    ~H5PartWriter() noexcept;

    //!  Write time stamp to ExodusII file
    void writeTimeStamp( tk::real time ) const;

  private:
    const std::string m_filename;          //!< File name
    H5PartFile* m_outFile;                 //!< H5Part file handle
};

} // tk::

#endif // H5PartWriter_h
