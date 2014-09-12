//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Thu 11 Sep 2014 04:30:54 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************
#ifndef PDF_h
#define PDF_h

#include <cmath>
#include <unordered_map>

#include <Types.h>
#include <PUPUtil.h>

namespace quinoa {

//! Univariate PDF estimator
class PDF {

    //! Univariate PDF is an unordered_map: key: bin id, value: sample counter
    using pdf = std::unordered_map< int, tk::real >;

  public:
    //! Empty constructor for Charm++
    explicit PDF() : m_binsize( 0 ), m_nsample( 0 ) {}

    //! Constructor
    //! \param[in]  binsize  Sample space bin size
    explicit PDF( tk::real binsize ) : m_binsize( binsize ), m_nsample( 0 ) {}

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    int nsample() const noexcept { return m_nsample; }

    //! Add new sample into univariate PDF
    //! \param[in]  sample  Value to insert
    void add( const tk::real& sample ) {
      ++m_nsample;
      ++m_pdf[ floor( sample / m_binsize + 0.5 ) ];
    }

    //! Add samples of a PDF passed in into PDF held
    //! \param[in]  p  New PDF to add
    void add( const PDF& p ) {
      m_binsize = p.binsize();
      m_nsample += p.nsample();
      for (const auto& e : p.map()) m_pdf[ e.first ] += e.second;
    }

    //! Zero bins
    void zero() noexcept { m_nsample = 0; m_pdf.clear(); }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf& map() const noexcept { return m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    tk::real binsize() const noexcept { return m_binsize; }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_nsample;
      using tk::operator|;
      p | m_pdf;
    }
    friend void operator|( PUP::er& p, PDF& t ) { t.pup(p); } 

  private:
    tk::real m_binsize;         //!< Sample space bin size
    int m_nsample;              //!< Number of samples collected
    pdf m_pdf;                  //!< Probability density function
};

} // quinoa::

#endif // PDF_h
