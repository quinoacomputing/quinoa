//******************************************************************************
/*!
  \file      src/Statistics/TriPDF.h
  \author    J. Bakosi
  \date      Thu 02 Oct 2014 12:26:02 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Joint trivariate PDF estimator
  \details   Joint trivariate PDF estimator
*/
//******************************************************************************
#ifndef TriPDF_h
#define TriPDF_h

#include <array>
#include <unordered_map>
#include <algorithm>

#include <Types.h>
#include <PUPUtil.h>

namespace quinoa {

//! Joint trivariate PDF estimator
class TriPDF {

  public:
    //! Number of sample space dimensions
    static const std::size_t dim = 3;

    //! Key type
    using key_type = std::array< long, dim >;

    //! Pair type
    using pair_type = std::pair< const key_type, tk::real >;

    // Hash function for key_type
    struct key_hash {
      long operator()( const key_type& key ) const {
        return std::hash< long >()( key[0] ) ^
               std::hash< long >()( key[1] ) ^
               std::hash< long >()( key[2] );
      }
    };

    //! Joint trivariate PDF is an unordered_map: key: three bin ids
    //! corresponding to the three sample space dimensions, mapped value: sample
    //! counter, hasher: XORed hash of the three bin ids
    using map_type = std::unordered_map< key_type, tk::real, key_hash >;

    //! Empty constructor for Charm++
    explicit TriPDF() : m_binsize( {{ 0, 0, 0 }} ), m_nsample( 0 ) {}

    //! Constructor: Initialize joint trivariate PDF container
    //! \param[in]  bs  Sample space bin size in all three directions
    explicit TriPDF( const std::vector< tk::real >& bs ) :
      m_binsize( {{ bs[0], bs[1], bs[2] }} ),
      m_nsample( 0 ) {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t nsample() const noexcept { return m_nsample; }

    //! Add sample to trivariate PDF
    //! \param[in]  sample  Sample to add
    void add( std::array< tk::real, dim > sample ) {
      ++m_nsample;
      ++m_pdf[ {{ std::lround( sample[0] / m_binsize[0] ),
                  std::lround( sample[1] / m_binsize[1] ),
                  std::lround( sample[2] / m_binsize[2] ) }} ];
    }

    //! Add multiple samples from a PDF
    //! \param[in]  p  PDF whose samples to add
    void addPDF( const TriPDF& p ) {
      m_binsize = p.binsize();
      m_nsample += p.nsample();
      for (const auto& e : p.map()) m_pdf[ e.first ] += e.second;
    }

    //! Zero bins
    void zero() noexcept { m_nsample = 0; m_pdf.clear(); }

    //! Constant accessor to underlying PDF map
    //! \return Reference to underlying map
    const map_type& map() const noexcept { return m_pdf; }

    //! Constant accessor to bin sizes
    //! \return Sample space bin sizes
    const std::array< tk::real, dim >& binsize() const noexcept
    { return m_binsize; }

    //! Return minimum and maximum bin ids of sample space in all three
    //! dimensions
    //! \return  {xmin,xmax,ymin,ymax,zmin,zmax}  Minima and maxima of bin ids
    std::array< long, 2*dim > extents() const {
      auto x = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first[0] < b.first[0]; } );
      auto y = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first[1] < b.first[1]; } );
      auto z = std::minmax_element( begin(m_pdf), end(m_pdf),
                 []( const pair_type& a, const pair_type& b )
                 { return a.first[2] < b.first[2]; } );
      return {{ x.first->first[0], x.second->first[0],
                y.first->first[1], y.second->first[1],
                z.first->first[2], z.second->first[2] }};
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      p | m_binsize;
      p | m_nsample;
      p | m_pdf;
    }

  private:
    std::array< tk::real, dim > m_binsize;   //!< Sample space bin sizes
    std::size_t m_nsample;                   //!< Number of samples collected
    map_type m_pdf;                          //!< Probability density function
};

} // quinoa::

#endif // TriPDF_h
