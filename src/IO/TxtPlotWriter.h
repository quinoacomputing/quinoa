//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.h
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 09:56:42 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Text plot writer
  \details   Text plot writer
*/
//******************************************************************************
#ifndef TxtPlotWriter_h
#define TxtPlotWriter_h

#include <string>

#include <QuinoaTypes.h>
#include <Writer.h>

namespace Quinoa {

class Statistics;

//! TxtPlotWriter : Writer
class TxtPlotWriter : public Writer {

  public:
    //! Constructor: Acquire plot file handle
    explicit TxtPlotWriter(const std::string& filename,
                           Statistics* const statistics);

    //! Destructor: Release plot file handle
    ~TxtPlotWriter() noexcept = default;

    //! Write out plot header
    void header();

    //! Write plot file
    void write(const int it, const real t);

  private:
    //! Don't permit copy constructor
    TxtPlotWriter(const TxtPlotWriter&) = delete;
    //! Don't permit copy assigment
    TxtPlotWriter& operator=(const TxtPlotWriter&) = delete;
    //! Don't permit move constructor
    TxtPlotWriter(TxtPlotWriter&&) = delete;
    //! Don't permit move assigment
    TxtPlotWriter& operator=(TxtPlotWriter&&) = delete;

    Statistics* const m_statistics;     //!< Statistics estimator
    const int m_nord;                   //!< Number of ordinary moments
    const int m_ncen;                   //!< Number of central moments
    const real* const m_ordinary;       //!< Ordinary moments
    const real* const m_central;        //!< Central moments
};

} // namespace Quinoa

#endif // TxtPlotWriter_h
