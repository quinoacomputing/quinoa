//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Tue May  7 12:17:57 2013
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

using namespace std;

namespace Quinoa {

//! Univariate PDF estimator
class PDF : public Distribution {

    //! Univariate PDF as unordered_map: key: bin id,
    //!                                  mapped value: sample counter
    using pdf = unordered_map<int,real>;

  public:
    //! Constructor: Initialize univariate PDF container
    explicit PDF(const real& binsize);

    //! Destructor: Clear univariate PDF container
    virtual ~PDF() noexcept;

    //! Insert new sample into univariate PDF
    virtual void insert(const real& sample);

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

} // namespace Quinoa

#endif // PDF_h
