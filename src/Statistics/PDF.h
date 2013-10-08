//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Mon 07 Oct 2013 08:47:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************
#ifndef PDF_h
#define PDF_h

#include <cmath>
#include <unordered_map>

#include <Types.h>
#include <Distribution.h>

namespace tk {

//! Univariate PDF estimator
class PDF : public Distribution {

    //! Univariate PDF as unordered_map: key: bin id,
    //!                                  mapped value: sample counter
    using pdf = std::unordered_map<int,tk::real>;

  public:
    //! Constructor: Initialize univariate PDF container
    //! \param[in]   binsize    Sample space bin size
    explicit PDF(const tk::real& binsize) : m_binsize(binsize), m_pdf() {}

    //! Destructor: Clear univariate PDF container
    ~PDF() noexcept override { m_pdf.clear(); }

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    const int& getNsample() const noexcept override { return m_nsample; }

    //! Insert new sample into univariate PDF
    //! \param[in]   sample    Value to insert
    void insert(const tk::real& sample) {
      ++m_nsample;                  // Increase number of samples in joint PDF
      ++m_pdf[floor(sample/m_binsize+0.5)];         // Add sample to joint PDF
    }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf* getMap() const noexcept { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const tk::real& getBinsize() const noexcept { return m_binsize; }

  private:
    //! Don't permit copy constructor
    PDF(const PDF&) = delete;
    //! Don't permit copy assigment
    PDF& operator=(const PDF&) = delete;
    //! Don't permit move constructor
    PDF(PDF&&) = delete;
    //! Don't permit move assigment
    PDF& operator=(PDF&&) = delete;

    const tk::real m_binsize;   //!< Sample space bin size
    pdf m_pdf;                  //!< Probability density function
};

} // tk::

#endif // PDF_h
