//******************************************************************************
/*!
  \file      src/IO/TxtPlotWriter.h
  \author    J. Bakosi
  \date      Tue May  7 13:08:36 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Text plot writer
  \details   Text plot writer
*/
//******************************************************************************
#ifndef TxtPlotWriter_h
#define TxtPlotWriter_h

#include <string>
#include <fstream>

#include <QuinoaTypes.h>
#include <PlotWriter.h>

using namespace std;

namespace Quinoa {

class Statistics;

//! TxtPlotWriter : PlotWriter
class TxtPlotWriter : public PlotWriter {

  public:
    //! Constructor: Acquire plot file handle
    explicit TxtPlotWriter(const string& filename,
                           Statistics* const statistics);

    //! Destructor: Release plot file handle
    ~TxtPlotWriter() noexcept = default;

    //! Write plot file
    virtual void write(const int it, const real t);

    //! Write out plot header
    void header();

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
