// *****************************************************************************
/*!
  \file      src/Statistics/BiPDF.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Joint bivariate PDF estimator
  \details   Joint bivariate PDF estimator. This class can be used to estimate a
    joint probability density function (PDF) of two scalar variables from an
    ensemble. The implementation uses the standard container std::unordered_map,
    which is a hash-based associative container with linear algorithmic
    complexity for insertion of a new sample.
*/
// *****************************************************************************
#ifndef BiPDF_h
#define BiPDF_h

#include <array>
#include <unordered_map>
#include <algorithm>

#include "Types.h"
#include "PUPUtil.h"

namespace tk {

//! Joint bivariate PDF estimator
class BiPDF {

  public:
    //! Number of sample space dimensions
    static const std::size_t dim = 2;

    //! Key type
    using key_type = std::array< long, dim >;

    //! Pair type
    using pair_type = std::pair< const key_type, tk::real >;

    // Hash functor for key_type
    struct key_hash {
      std::size_t operator()( const key_type& key ) const {
        return std::hash< long >()( key[0] ) ^ std::hash< long >()( key[1] );
      }
    };

    //! \brief Joint bivariate PDF
    //! \details The underlying container type is an unordered_map where the key
    //!   is two bin ids corresponding to the two sample space dimensions, and
    //!   the mapped value is the sample counter. The hasher functor, defined by
    //!   key_hash provides an XORed hash of the two bin ids.
    using map_type = std::unordered_map< key_type, tk::real, key_hash >;

    //! Empty constructor for Charm++
    explicit BiPDF() : m_binsize( {{ 0, 0 }} ), m_nsample( 0 ), m_pdf() {}

    //! Constructor: Initialize joint bivariate PDF container
    //! \param[in] bs Sample space bin size in both directions
    explicit BiPDF( const std::vector< tk::real >& bs ) :
      m_binsize( {{ bs[0], bs[1] }} ), m_nsample( 0 ), m_pdf() {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t nsample() const noexcept { return m_nsample; }

    //! Add sample to bivariate PDF
    //! \param[in] sample Sample to add
    void add( std::array< tk::real, dim > sample ) {
      ++m_nsample;
      ++m_pdf[ {{ std::lround( sample[0] / m_binsize[0] ),
                  std::lround( sample[1] / m_binsize[1] ) }} ];
    }

    //! Add multiple samples from a PDF
    //! \param[in] p PDF whose samples to add
    void addPDF( const BiPDF& p ) {
      m_binsize = p.binsize();
      m_nsample += p.nsample();
      for (const auto& e : p.map()) m_pdf[ e.first ] += e.second;
    }

    //! Zero bins
    void zero() noexcept { m_nsample = 0; m_pdf.clear(); }

    //! Constant accessor to underlying PDF map
    //! \return Constant reference to underlying map
    const map_type& map() const noexcept { return m_pdf; }

    //! Constant accessor to bin sizes
    //! \return Constant reference to sample space bin sizes
    const std::array< tk::real, dim >& binsize() const noexcept
    { return m_binsize; }

    //! Return minimum and maximum bin ids of sample space in both dimensions
    //! \return {xmin,xmax,ymin,ymax} Minima and maxima of the bin ids in a
    //!    std::array
    std::array< long, 2*dim > extents() const {
      auto x = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first[0] < b.first[0]; } );
      auto y = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first[1] < b.first[1]; } );
      return {{ x.first->first[0], x.second->first[0],
                y.first->first[1], y.second->first[1] }};
    }

    /** @name Pack/Unpack: Serialize BiPDF object for Charm++ */
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
    //! \param[in,out] c BiPDF object reference
    friend void operator|( PUP::er& p, BiPDF& c ) { c.pup(p); }
    ///@}

  private:
    std::array< tk::real, dim > m_binsize;  //!< Sample space bin sizes
    std::size_t m_nsample;                  //!< Number of samples collected
    map_type m_pdf;                         //!< Probability density function
};

} // tk::

#endif // BiPDF_h
