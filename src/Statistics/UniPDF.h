//******************************************************************************
/*!
  \file      src/Statistics/UniPDF.h
  \author    J. Bakosi
  \date      Sat 27 Sep 2014 09:54:54 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator
*/
//******************************************************************************
#ifndef UniPDF_h
#define UniPDF_h

#include <unordered_map>
#include <algorithm>

#include <Types.h>
#include <PUPUtil.h>

namespace quinoa {

//! Univariate PDF estimator
class UniPDF {

  public:
    //! Key type
    using key_type = long;

    //! Pair type
    using pair_type = std::pair< const key_type, tk::real >;

    //! Univariate PDF is an unordered_map: key: bin id, value: sample counter
    using map_type = std::unordered_map< key_type, tk::real >;

    //! Empty constructor for Charm++
    explicit UniPDF() : m_binsize( 0 ), m_nsample( 0 ) {}

    //! Constructor
    //! \param[in]  binsize  Sample space bin size
    explicit UniPDF( tk::real binsize ) :
      m_binsize( binsize ), m_nsample( 0 ) {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t nsample() const noexcept { return m_nsample; }

    //! Add sample to univariate PDF
    //! \param[in]  sample  Sample to insert
    void add( tk::real sample ) {
      ++m_nsample;
      ++m_pdf[ std::lround( sample / m_binsize ) ];
    }

    //! Add multiple samples from a PDF
    //! \param[in]  p  PDF whose samples to add
    void addPDF( const UniPDF& p ) {
      m_binsize = p.binsize();
      m_nsample += p.nsample();
      for (const auto& e : p.map()) m_pdf[ e.first ] += e.second;
    }

    //! Zero bins
    void zero() noexcept { m_nsample = 0; m_pdf.clear(); }

    //! Constant accessor to underlying PDF map
    //! \return Constant reference to underlying map
    const map_type& map() const noexcept { return m_pdf; }

    //! Constant accessor to bin size
    //! \return Sample space bin size
    tk::real binsize() const noexcept { return m_binsize; }

    //! Return minimum and maximum bin ids of sample space
    //! \return  {min,max}  Minimum and maximum of the bin ids
    std::pair< key_type, key_type > extents() const {
      auto x = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first < b.first; } );
      return { { x.first->first }, { x.second->first } };
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_nsample;
      p | m_pdf;
    }

  private:
    tk::real m_binsize;         //!< Sample space bin size
    std::size_t m_nsample;      //!< Number of samples collected
    map_type m_pdf;             //!< Probability density function
};

} // quinoa::

#endif // UniPDF_h
