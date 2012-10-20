//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 09:49:03 PM MDT
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

//! PDF estimator base class
class PDF {

  typedef unordered_map<int,real> Pdf;

  public:
    //! Constructor: Initialize PDF container
    PDF(const real binsize) : m_binsize(binsize), m_sample(0) {}

    //! Destructor: Clear PDF container
    ~PDF() { m_pdf.clear(); }

    //! Insert new value into PDF
    //! \param[in]   value    Value to insert
    void insert(const real& value) {
      ++m_sample;
      ++m_pdf[floor(value/m_binsize+0.5)];
    }

    //! Save PDF to file
    void save(const string& filename);

  private:
    //! Don't permit copy constructor
    PDF(const PDF&) = delete;
    //! Don't permit copy assigment
    PDF& operator=(const PDF&) = delete;
    //! Don't permit move constructor
    PDF(PDF&&) = delete;
    //! Don't permit move assigment
    PDF& operator=(PDF&&) = delete;

    const real m_binsize;   //!< Bin size = (max-min)/nbin
    int m_sample;           //!< Number of samples
    Pdf m_pdf;              //!< Probability density function
};

} // namespace Quinoa

#endif // PDF_h
