//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 10:59:14 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF estimator base class
  \details   PDF estimator base class
*/
//******************************************************************************
#ifndef PDF_h
#define PDF_h

#include <cmath>
#include <unordered_map>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! PDF as unordered_map: bin id as key, counter as mapped value
typedef unordered_map<int,real> Pdf;

//! PDF estimator base class
class PDF {

  public:
    //! Constructor: Initialize PDF container
    //! \param[in]   binsize    Sample space bin size
    PDF(const real binsize) : m_binsize(binsize), m_nsample(0) {}

    //! Destructor: Clear PDF container
    ~PDF() { m_pdf.clear(); }

    //! Insert new value into PDF
    //! \param[in]   value    Value to insert
    void insert(const real& value) {
      ++m_nsample;
      ++m_pdf[floor(value/m_binsize+0.5)];
    }

    //! Constant accessor to PDF
    //! \param[out] Pointer to Pdf
    const Pdf* getPDF() const { return &m_pdf; }
    //! Constant accessor to binsize
    //! \param[out] Sample space bin size
    const real& getBinsize() const { return m_binsize; }
    //! Constant accessor to number of samples
    //! \param[out] Number of samples collected
    const int& getNsample() const { return m_nsample; }

  private:
    //! Don't permit copy constructor
    PDF(const PDF&) = delete;
    //! Don't permit copy assigment
    PDF& operator=(const PDF&) = delete;
    //! Don't permit move constructor
    PDF(PDF&&) = delete;
    //! Don't permit move assigment
    PDF& operator=(PDF&&) = delete;

    const real m_binsize;   //!< Sample space bin size
    int m_nsample;          //!< Number of samples collected
    Pdf m_pdf;              //!< Probability density function
};

} // namespace Quinoa

#endif // PDF_h
