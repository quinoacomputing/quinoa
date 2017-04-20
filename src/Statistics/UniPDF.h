// *****************************************************************************
/*!
  \file      src/Statistics/UniPDF.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator. This class can be used to estimate a
    probability density function of (PDF) a scalar variable from an ensemble.
    The implementation uses the standard container std::unordered_map, which is
    a hash-based associative container with linear algorithmic complexity for
    insertion of a new sample.
*/
// *****************************************************************************
#ifndef UniPDF_h
#define UniPDF_h

#include <array>
#include <unordered_map>
#include <algorithm>

#include "Types.h"
#include "Exception.h"
#include "PUPUtil.h"

namespace tk {

//! Univariate PDF estimator
class UniPDF {

  public:
    //! Number of sample space dimensions
    static const std::size_t dim = 1;

    //! Key type
    using key_type = long;

    //! Pair type
    using pair_type = std::pair< const key_type, tk::real >;

    //! \brief Univariate PDF
    //! \details The underlying container type is an unordered_map where the key
    //!   is one bin id corresponding to the single sample space dimension, and
    //!   the mapped value is the sample counter. The hasher functor used here
    //!   is the default for the key type provided by the standard library.
    using map_type = std::unordered_map< key_type, tk::real >;

    //! Empty constructor for Charm++
    explicit UniPDF() : m_binsize( 0 ), m_nsample( 0 ), m_pdf() {}

    //! Constructor: Initialize univariate PDF container
    //! \param[in] bs Sample space bin size
    explicit UniPDF( tk::real bs ) :
      m_binsize( bs ), m_nsample( 0 ), m_pdf() {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t nsample() const noexcept { return m_nsample; }

    //! Add sample to univariate PDF
    //! \param[in] sample Sample to insert
    void add( tk::real sample ) {
      Assert( m_binsize > 0, "Bin size must be positive" );
      ++m_nsample;
      ++m_pdf[ std::lround( sample / m_binsize ) ];
    }

    //! Add multiple samples from a PDF
    //! \param[in] p PDF whose samples to add
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
    //! \return {min,max} Minimum and maximum of the bin ids
    std::array< long, 2*dim > extents() const {
      auto x = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first < b.first; } );
      return {{ x.first->first, x.second->first }};
    }

    /** @name Pack/Unpack: Serialize UniPDF object for Charm++ */
    ///@{
    //! Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_nsample;
      p | m_pdf;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c UniPDF object reference
    friend void operator|( PUP::er& p, UniPDF& c ) { c.pup(p); }
    ///@}

  private:
    tk::real m_binsize;         //!< Sample space bin size
    std::size_t m_nsample;      //!< Number of samples collected
    map_type m_pdf;             //!< Probability density function
};

//! Output univariate PDF to output stream
//! \param[in,out] os Stream to output to
//! \param[in] p PDF to output
//! \return Updated stream
//! \note Used for debugging.
//! \author J. Bakosi
static inline
std::ostream& operator<< ( std::ostream& os, const tk::UniPDF& p ) {
  os << p.binsize() << ", " << p.nsample() << ": ";
  std::map< typename tk::UniPDF::key_type, tk::real >
    sorted( p.map().begin(), p.map().end() );
  for (const auto& b : sorted) os << '(' << b.first << ',' << b.second << ") ";
  return os;
}

} // tk::

#endif // UniPDF_h
