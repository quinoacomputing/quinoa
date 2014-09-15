//******************************************************************************
/*!
  \file      src/Statistics/JPDF.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 08:38:58 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************
#ifndef JPDF_h
#define JPDF_h

#include <vector>
#include <unordered_map>
#include <algorithm>

#include <Types.h>

namespace quinoa {

//! Joint PDF estimator
class JPDF {

  public:
    //! Key type
    using key_type = std::vector< std::size_t >;

    // Hash function for std::vector<int>
    struct key_hash {
      size_t operator()(const key_type& key) const {
        size_t h = 0;
        for (auto& k : key) h ^= std::hash< std::size_t >()( k );
        return h;
      }
    };

    //! Joint PDF as unordered_map: key: bin ids,
    //                              mapped value: sample counter,
    //                              hasher: XORed hash of all bin ids
    using pdf = std::unordered_map< key_type, tk::real, key_hash >;

    //! Constructor: Initialize joint PDF container
    //! \param[in]   dim        Dimension of sample space
    //! \param[in]   binsize    Sample space bin size
    explicit JPDF( int dim, tk::real binsize ) :
      m_binsize( binsize ),
      m_size( 0 ),
      m_key( dim ) {}

    //! Accessor to number of samples
    //! \return Number of samples collected
    std::size_t size() const noexcept { return m_size; }

    //! Insert new sample into joint PDF
    void insert(const std::vector< tk::real >& sample);

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf& map() const noexcept { return m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    tk::real binsize() const noexcept { return m_binsize; }

  private:
    tk::real m_binsize;      //!< Sample space bin size
    std::size_t m_size;      //!< Number of samples collected
    key_type m_key;          //!< Temporary key for finding the sample space bin
    pdf m_pdf;               //!< Probability density function
};

} // quinoa::

#endif // JPDF_h
