//******************************************************************************
/*!
  \file      src/IO/TxtStatWriter.h
  \author    J. Bakosi
  \date      Sat 24 Aug 2013 08:47:02 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Text statistics writer
  \details   Text statistics writer
*/
//******************************************************************************
#ifndef TxtStatWriter_h
#define TxtStatWriter_h

#include <string>

#include <QuinoaTypes.h>
#include <Writer.h>

namespace Quinoa {

class Statistics;

//! TxtStatWriter : Writer
class TxtStatWriter : public Writer {

  public:
    //! Constructor
    explicit TxtStatWriter(const std::string& filename,
                           Statistics* const statistics);

    //! Destructor
    virtual ~TxtStatWriter() noexcept = default;

    //! Write out statistics file header
    void header();

    //! Write statistics file
    void write(const int it, const real t);

  private:
    //! Don't permit copy constructor
    TxtStatWriter(const TxtStatWriter&) = delete;
    //! Don't permit copy assigment
    TxtStatWriter& operator=(const TxtStatWriter&) = delete;
    //! Don't permit move constructor
    TxtStatWriter(TxtStatWriter&&) = delete;
    //! Don't permit move assigment
    TxtStatWriter& operator=(TxtStatWriter&&) = delete;

    Statistics* const m_statistics;     //!< Statistics estimator
    const int m_nord;                   //!< Number of ordinary moments
    const int m_ncen;                   //!< Number of central moments
    const real* const m_ordinary;       //!< Ordinary moments
    const real* const m_central;        //!< Central moments
};

} // namespace Quinoa

#endif // TxtStatWriter_h
