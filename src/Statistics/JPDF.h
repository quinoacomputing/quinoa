//******************************************************************************
/*!
  \file      src/Statistics/JPDF.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 08:00:56 PM MST
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
#include <Distribution.h>
#include <StatException.h>

using namespace std;

namespace Quinoa {

//! Joint PDF estimator
class JPDF : private Distribution {

  private:
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

  public:
    //! Constructor: Initialize joint PDF container
    JPDF(const int dim, const real binsize);

    //! Destructor: Clear joint PDF container
    virtual ~JPDF();

    //! Insert new value into joint PDF
    virtual void insert(const vector<real>& value);

    //! Constant accessor to number of samples
    //! \return Number of samples collected
    const int& getNsample() const { return m_nsample; }

    //! Constant accessor to PDF map
    //! \return Pointer to map
    const Pdf* getMap() const { return &m_pdf; }

    //! Constant accessor to binsize
    //! \return Sample space bin size
    const real& getBinsize() const { return m_binsize; }

  private:
    //! Don't permit copy constructor
    JPDF(const JPDF&) = delete;
    //! Don't permit copy assigment
    JPDF& operator=(const JPDF&) = delete;
    //! Don't permit move constructor
    JPDF(JPDF&&) = delete;
    //! Don't permit move assigment
    JPDF& operator=(JPDF&&) = delete;

    key_type m_key;         //!< Temporary key for finding the sample space bin
    const int m_dim;        //!< Sample space dimension
    const real m_binsize;   //!< Sample space bin size
    Pdf m_pdf;              //!< Probability density function
};

} // namespace Quinoa

#endif // JPDF_h
