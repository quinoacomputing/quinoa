// *****************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:26:15 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "Types.h"
#include "Writer.h"
#include "Keywords.h"
#include "Options/TxtFloatFormat.h"

namespace tk {

//! \brief TxtStatWriter : tk::Writer
//! \details ASCII statistics writer class that facilitates outputing
//!   statistics to text files.
//! \author J. Bakosi
class TxtStatWriter : public tk::Writer {

  public:
    //! Constructor
    explicit TxtStatWriter(
      const std::string& filename,
      tk::ctr::TxtFloatFormatType format = tk::ctr::TxtFloatFormatType::DEFAULT,
      kw::precision::info::expect::type precision = std::cout.precision(),
      std::ios_base::openmode mode = std::ios_base::out );

    //! Write out statistics file header
    void header( const std::vector< std::string >& nameOrd,
                 const std::vector< std::string >& nameCen ) const;

    //! Write statistics file
    std::size_t stat( uint64_t it,
                      tk::real t,
                      const std::vector< tk::real >& ordinary,
                      const std::vector< tk::real >& central );

  private:
    int m_precision;    //!< Floating-point precision in digits
    int m_width;        //!< Floating-point number width
};

} // tk::

#endif // TxtStatWriter_h
