//******************************************************************************
/*!
  \file      src/Statistics/PDF.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 08:42:26 AM MDT
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
    using pdf = std::unordered_map< std::size_t, tk::real >;

  public:
    //! Empty constructor for Charm++
    explicit PDF() : m_binsize( 0 ), m_size( 0 ) {}

    //! Constructor
    //! \param[in]  binsize  Sample space bin size
    explicit PDF( tk::real binsize ) : m_binsize( binsize ), m_size( 0 ) {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t size() const noexcept { return m_size; }

    //! Add new sample into univariate PDF
    //! \param[in]  sample  Value to insert
    void add( tk::real sample ) {
      ++m_size;
      ++m_pdf[ floor( sample / m_binsize + 0.5 ) ];
    }

    //! Add samples of a PDF passed in into PDF held
    //! \param[in]  p  New PDF to add
    void add( const PDF& p ) {
      m_binsize = p.binsize();
      m_size += p.size();
      for (const auto& e : p.map()) m_pdf[ e.first ] += e.second;
    }

    //! Zero bins
    void zero() noexcept { m_size = 0; m_pdf.clear(); }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf& map() const noexcept { return m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    tk::real binsize() const noexcept { return m_binsize; }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_size;
      tk::pup( p, m_pdf );
    }

  private:
    tk::real m_binsize;         //!< Sample space bin size
    std::size_t m_size;         //!< Number of samples collected
    pdf m_pdf;                  //!< Probability density function
};

} // quinoa::

#endif // PDF_h
