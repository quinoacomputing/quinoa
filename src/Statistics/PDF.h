//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Sat 27 Oct 2012 11:45:41 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************
#ifndef PDF_h
#define PDF_h

#include <cmath>
#include <unordered_map>

#include <QuinoaTypes.h>
#include <Distribution.h>
#include <StatException.h>

using namespace std;

namespace Quinoa {

//! Univariate PDF estimator
class PDF : Distribution {

    //! Univariate PDF as unordered_map: key: bin id,
    //!                                  mapped value: sample counter
    using Pdf = unordered_map<int,real>;

  public:
    //! Constructor: Initialize univariate PDF container
    //! \param[in]   binsize    Sample space bin size
    PDF(const real& binsize) : m_binsize(binsize) {}

    //! Destructor: Clear univariate PDF container
    ~PDF() { m_pdf.clear(); }

    //! Insert new value into univariate PDF
    //! \param[in]   value    Value to insert
    void insert(const real& value) {
      ++m_nsample;
      ++m_pdf[floor(value/m_binsize+0.5)];
    }
    //! Throw exception if vector sample is given
    void insert(const vector<real>& value) {
      throw StatException(WARNING, STATEXCEPT_UNIMPLEMENTED);
    }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const Pdf* getMap() const { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const real& getBinsize() const { return m_binsize; }

    //! Constant accessor to number of samples
    //! \return Number of samples collected
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
    Pdf m_pdf;              //!< Probability density function
};

} // namespace Quinoa

#endif // PDF_h
