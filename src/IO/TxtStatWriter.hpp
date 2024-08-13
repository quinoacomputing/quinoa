// *****************************************************************************
/*!
  \file      src/IO/TxtStatWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Text statistics writer declaration
  \details   This file declares the ASCII statistics writer class that
     facilitates outputing statistics to text files.
*/
// *****************************************************************************
#ifndef TxtStatWriter_h
#define TxtStatWriter_h

#include <string>
#include <vector>
#include <fstream>

#include "Types.hpp"
#include "Writer.hpp"
#include "Options/TxtFloatFormat.hpp"

namespace tk {

//! \brief TxtStatWriter : tk::Writer
//! \details ASCII statistics writer class that facilitates outputing
//!   statistics to text files.
class TxtStatWriter : public tk::Writer {

  public:
    //! Constructor
    explicit TxtStatWriter(
      const std::string& filename,
      tk::ctr::TxtFloatFormatType format = tk::ctr::TxtFloatFormatType::DEFAULT,
      std::streamsize precision = std::cout.precision(),
      std::ios_base::openmode mode = std::ios_base::out );

    //! Write out statistics file header
    void header( const std::vector< std::string >& nameOrd,
                 const std::vector< std::string >& nameCen,
                 const std::vector< std::string >& nameExt ) const;

    //! Write statistics file
    std::size_t stat( uint64_t it,
                      tk::real t,
                      const std::vector< tk::real >& ordinary,
                      const std::vector< tk::real >& central,
                      const std::vector< tk::real >& extra );

  private:
    int m_precision;    //!< Floating-point precision in digits
    int m_width;        //!< Floating-point number width
};

} // tk::

#endif // TxtStatWriter_h
