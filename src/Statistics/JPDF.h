//******************************************************************************
/*!
  \file      src/Statistics/JPDF.h
  \author    J. Bakosi
  \date      Mon 07 Oct 2013 08:46:12 PM MDT
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
#include <Distribution.h>

namespace tk {

//! Joint PDF estimator
class JPDF : public Distribution {

  public:
    //! Key type
    using key_type = std::vector<int>;

    // Hash function for std::vector<int>
    struct key_hash {
      size_t operator()(const key_type& key) const {
        size_t h = 0;
        for (auto& k : key) h ^= std::hash<int>()(k);
        return h;
      }
    };

    //! Joint PDF as unordered_map: key: bin ids,
    //                              mapped value: sample counter,
    //                              hasher: XORed hash of all bin ids
    using pdf = std::unordered_map<key_type, tk::real, key_hash>;

    //! Constructor: Initialize joint PDF container
    //! \param[in]   dim        Dimension of sample space
    //! \param[in]   binsize    Sample space bin size
    explicit JPDF(const int dim, const tk::real binsize) :
      m_binsize(binsize),
      m_key(dim),
      m_pdf() {}

    //! Destructor: Clear joint PDF container
    ~JPDF() noexcept override { m_pdf.clear(); }

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    const int& getNsample() const noexcept override { return m_nsample; }

    //! Insert new sample into joint PDF
    void insert(const std::vector<tk::real>& sample);

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf* getMap() const noexcept { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const tk::real& getBinsize() const noexcept { return m_binsize; }

  private:
    //! Don't permit copy constructor
    JPDF(const JPDF&) = delete;
    //! Don't permit copy assigment
    JPDF& operator=(const JPDF&) = delete;
    //! Don't permit move constructor
    JPDF(JPDF&&) = delete;
    //! Don't permit move assigment
    JPDF& operator=(JPDF&&) = delete;

    const tk::real m_binsize;//!< Sample space bin size
    key_type m_key;          //!< Temporary key for finding the sample space bin
    pdf m_pdf;               //!< Probability density function
};

} // tk::

#endif // JPDF_h
