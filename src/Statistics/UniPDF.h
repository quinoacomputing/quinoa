//******************************************************************************
/*!
  \file      src/Statistics/UniPDF.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:08:28 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Univariate PDF estimator
  \details   Univariate PDF estimator. This class can be used to estimate a
    probability density function of (PDF) a scalar variable from an ensemble.
    The implementation uses the standard container std::unordered_map, which is
    a hash-based associative container with linear algorithmic complexity for
    insertion of a new sample.
*/
//******************************************************************************
#ifndef UniPDF_h
#define UniPDF_h

#include <array>
#include <unordered_map>
#include <algorithm>

#include "Types.h"
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
    explicit UniPDF() : m_binsize( 0 ), m_nsample( 0 ) {}

    //! Constructor: Initialize univariate PDF container
    //! \param[in] binsize Sample space bin size
    explicit UniPDF( tk::real binsize ) :
      m_binsize( binsize ), m_nsample( 0 ) {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t nsample() const noexcept { return m_nsample; }

    //! Add sample to univariate PDF
    //! \param[in] sample Sample to insert
    void add( tk::real sample ) {
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

    /** @name Pack/Unpack: Serialize BiPDF object for Charm++ */
    ///@{
    //! Pack/Unpack serialize member function
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_nsample;
      p | m_pdf;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    //! \param[inout] c UniPDF object reference
    friend void operator|( PUP::er& p, UniPDF& c ) { c.pup(p); }
    ///@}

  private:
    tk::real m_binsize;         //!< Sample space bin size
    std::size_t m_nsample;      //!< Number of samples collected
    map_type m_pdf;             //!< Probability density function
};

} // tk::

#endif // UniPDF_h
