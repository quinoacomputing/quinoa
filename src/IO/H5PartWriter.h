// *****************************************************************************
/*!
  \file      src/IO/H5PartWriter.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     H5Part particles data writer
  \details   H5Part particles data writer class declaration, facilitating
    writing particle coordinates and associated particle fields into HDF5-based
    H5Part data files in parallel, using MPI-IO.
*/
// *****************************************************************************
#ifndef H5PartWriter_h
#define H5PartWriter_h

#include <vector>
#include <string>
#include <memory>

#include "NoWarning/H5Part.h"

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

    //! Write particle coordinates to H5Part file
    void writeCoords( uint64_t it,
                      const std::vector< tk::real >& x,
                      const std::vector< tk::real >& y,
                      const std::vector< tk::real >& z ) const;

  private:
    const std::string m_filename;               //!< File name
};

} // tk::

#endif // H5PartWriter_h
