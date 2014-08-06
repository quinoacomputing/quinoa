//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 04:46:43 PM MDT
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
    explicit TxtStatWriter( const std::string& filename,
                            int nord, int ncen,
                            tk::real* ordinary, tk::real* central,
                            const std::vector< bool >& plotOrdinary,
                            const std::vector< std::string >& nameOrdinary,
                            const std::vector< std::string >& nameCentral ) :
      Writer( filename ),
      m_nord( nord ),
      m_ncen( ncen ),
      m_ordinary( ordinary ),
      m_central( central ),
      m_plotOrdinary( plotOrdinary ),
      m_nameOrdinary( nameOrdinary ),
      m_nameCentral( nameCentral ) {}

    //! Destructor
    ~TxtStatWriter() noexcept override = default;

    //! Write out statistics file header
    void header() const;

    //! Write statistics file
    void writeStat(const int it, const tk::real t);

  private:
    //! Don't permit copy constructor
    TxtStatWriter(const TxtStatWriter&) = delete;
    //! Don't permit copy assigment
    TxtStatWriter& operator=(const TxtStatWriter&) = delete;
    //! Don't permit move constructor
    TxtStatWriter(TxtStatWriter&&) = delete;
    //! Don't permit move assigment
    TxtStatWriter& operator=(TxtStatWriter&&) = delete;

    const int m_nord;                   //!< Number of ordinary moments
    const int m_ncen;                   //!< Number of central moments
    const tk::real* const m_ordinary;   //!< Ordinary moments
    const tk::real* const m_central;    //!< Central moments
    const std::vector< bool >& m_plotOrdinary;  //!< Whether to plot ord moments
    const std::vector< std::string >& m_nameOrdinary; //!< Ordinary moment names
    const std::vector< std::string >& m_nameCentral;  //!< Central moment names
};

} // quinoa::

#endif // TxtStatWriter_h
