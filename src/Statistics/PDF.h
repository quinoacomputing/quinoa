//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Wed 17 Oct 2012 09:18:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PDF estimator base class
  \details   PDF estimator base class
*/
//******************************************************************************
#ifndef PDF_h
#define PDF_h

#include <unordered_map>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! PDF estimator base class
class PDF {

  typedef unordered_map<Real,Real> Pdf;

  protected:
    //! Constructor: Initialize PDF container
    PDF(const Real min, const Real max, const Real nbin);

    //! Destructor: Clear PDF container
    ~PDF();

  private:
    //! Don't permit copy constructor
    PDF(const PDF&) = delete;
    //! Don't permit copy assigment
    PDF& operator=(const PDF&) = delete;
    //! Don't permit move constructor
    PDF(PDF&&) = delete;
    //! Don't permit move assigment
    PDF& operator=(PDF&&) = delete;

    const Real m_min;         //!< Most negative value of sample space
    const Real m_max;         //!< Most positive value of sample space
    const Real m_nbin;        //!< Number of bins of the sample space
    const Real m_binSize;     //!< Bin size = (max-min)/nbin
    const Real m_halfBinSize; //!< Bin size = (max-min)/nbin

    Pdf m_pdf;             //!< Probability density function
};

} // namespace Quinoa

#endif // PDF_h
