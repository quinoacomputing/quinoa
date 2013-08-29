//******************************************************************************
/*!
  \file      src/Statistics/JPDF.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:26:04 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************
#ifndef JPDF_h
#define JPDF_h

#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>

#include <QuinoaTypes.h>
#include <Distribution.h>

namespace quinoa {

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
    using pdf = std::unordered_map<key_type, real, key_hash>;
    //! Ordered counterpart
    using ordered_pdf = std::map<key_type, real, key_hash>;

    //! Constructor: Initialize joint PDF container
    //! \param[in]   dim        Dimension of sample space
    //! \param[in]   binsize    Sample space bin size
    explicit JPDF(const int dim, const real binsize) :
      m_binsize(binsize),
      m_key(dim),
      m_pdf() {}

    //! Destructor: Clear joint PDF container
    virtual ~JPDF() noexcept { m_pdf.clear(); }

    //! Insert new sample into joint PDF
    virtual void insert(const std::vector<real>& sample);

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    const int& getNsample() const noexcept { return m_nsample; }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const pdf* getMap() const noexcept { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const real& getBinsize() const noexcept { return m_binsize; }

  private:
    //! Don't permit copy constructor
    JPDF(const JPDF&) = delete;
    //! Don't permit copy assigment
    JPDF& operator=(const JPDF&) = delete;
    //! Don't permit move constructor
    JPDF(JPDF&&) = delete;
    //! Don't permit move assigment
    JPDF& operator=(JPDF&&) = delete;

    const real m_binsize;   //!< Sample space bin size
    key_type m_key;         //!< Temporary key for finding the sample space bin
    pdf m_pdf;              //!< Probability density function
};

} // namespace quinoa

#endif // JPDF_h
