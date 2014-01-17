//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 10:25:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRNG_h
#define MKLRNG_h

#include <memory>

#include <mkl_vsl_types.h>

#include <Options/RNG.h>
#include <RNG.h>

namespace tk {

//! MKL-based random number generator
class MKLRNG : public RNG {

  public:
    //! Constructor
    explicit MKLRNG( int nthreads,
                     int brng,
                     unsigned int seed,
                     int uniform_method,
                     int gaussian_method );

    //! Destructor
    ~MKLRNG() noexcept override;

    //! Uniform RNG
    void uniform(int tid, int num, double* r) const override;

    //! Gaussian RNG
    void gaussian(int tid, int num, double* r) const override;

  private:
    //! Don't permit copy constructor
    MKLRNG(const MKLRNG&) = delete;
    //! Don't permit copy assigment
    MKLRNG& operator=(const MKLRNG&) = delete;
    //! Don't permit move constructor
    MKLRNG(MKLRNG&&) = delete;
    //! Don't permit move assigment
    MKLRNG& operator=(MKLRNG&&) = delete;

    const int m_uniform_method;         //!< Uniform method to use
    const int m_gaussian_method;        //!< Gaussian method to use
    const int m_nthreads;               //!< Number of threads

    //! Random number stream for threads
    std::unique_ptr< VSLStreamStatePtr[] > m_stream;
};

} // tk::

#endif // MKLRNG_h
