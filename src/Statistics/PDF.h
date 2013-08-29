//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:26:11 2013
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

namespace quinoa {

//! Univariate PDF estimator
class PDF : public Distribution {

    //! Univariate PDF as unordered_map: key: bin id,
    //!                                  mapped value: sample counter
    using pdf = std::unordered_map<int,real>;

  public:
    //! Constructor: Initialize univariate PDF container
    //! \param[in]   binsize    Sample space bin size
    explicit PDF(const real& binsize) : m_binsize(binsize), m_pdf() {}

    //! Destructor: Clear univariate PDF container
    virtual ~PDF() noexcept { m_pdf.clear(); }

    //! Insert new sample into univariate PDF
    //! \param[in]   sample    Value to insert
    virtual void insert(const real& sample) {
      ++m_nsample;                  // Increase number of samples in joint PDF
      ++m_pdf[floor(sample/m_binsize+0.5)];         // Add sample to joint PDF
    }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf* getMap() const noexcept { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const real& getBinsize() const noexcept { return m_binsize; }

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    const int& getNsample() const noexcept { return m_nsample; }

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
    pdf m_pdf;              //!< Probability density function
};

} // namespace quinoa

#endif // PDF_h
