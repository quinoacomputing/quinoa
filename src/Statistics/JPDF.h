//******************************************************************************
/*!
  \file      src/Statistics/JPDF.h
  \author    J. Bakosi
  \date      Thu 25 Oct 2012 06:34:48 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Joint PDF estimator
  \details   Joint PDF estimator
*/
//******************************************************************************
#ifndef JPDF_h
#define JPDF_h

#include <vector>
#include <unordered_map>
#include <algorithm>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! Joint PDF estimator
class JPDF {

  public:
    //! Constructor: Initialize joint PDF container
    JPDF(const real dim, const real binsize);

    //! Destructor: Clear joint PDF container
    ~JPDF();

    //! Insert new value into joint PDF
    void insert(const vector<real>& value);

    //! Constant accessor to PDF
    //! \param[out] Pointer to Pdf
    //const Pdf* getPDF() const { return &m_pdf; }

    //! Constant accessor to binsize
    //! \param[out] Sample space bin size
    const real& getBinsize() const { return m_binsize; }

    //! Constant accessor to number of samples
    //! \param[out] Number of samples collected
    const int& getNsample() const { return m_nsample; }

  private:
    //! Don't permit copy constructor
    JPDF(const JPDF&) = delete;
    //! Don't permit copy assigment
    JPDF& operator=(const JPDF&) = delete;
    //! Don't permit move constructor
    JPDF(JPDF&&) = delete;
    //! Don't permit move assigment
    JPDF& operator=(JPDF&&) = delete;

    //! Key type
    using key_type = vector<int>;

    // Hash function for vector<int>
    struct key_hash {
      size_t operator()(const key_type& key) const {
        size_t h = 0;
        for (auto& k : key) h ^= hash<int>()(k);
        return h;
      }
    };

    //! Joint PDF as unordered_map: key: bin ids,
    //                              mapped value: sample counter,
    //                              hasher: XORed hash of all bin ids
    using Pdf = unordered_map<key_type,real,key_hash>;

    key_type m_key;         //!< Temporary key for finding the sample space bin
    const int m_dim;        //!< Sample space dimension
    const real m_binsize;   //!< Sample space bin size
    int m_nsample;          //!< Number of samples collected
    Pdf m_pdf;              //!< Probability density function
};

} // namespace Quinoa

#endif // JPDF_h
