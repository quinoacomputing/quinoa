//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.h
  \author    J. Bakosi
  \date      Fri 25 Oct 2013 10:36:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRNG_h
#define MKLRNG_h

#include <memory>

#include <mkl_vsl_types.h>

#include <RNG.h>

namespace quinoa {

//! MKL-based random number generator
class MKLRNG : public tk::RNG {

  public:
    //! Constructor
    explicit MKLRNG(const int nthreads, int brng, unsigned int seed);

    //! Destructor: Free all random number tables and streams
    virtual ~MKLRNG() noexcept;

    //! Uniform RNG
    void uniform(int tid, int num, tk::real* r) const override;

    //! Gaussian RNG
    void gaussian(int tid, int num, tk::real* r) const override;

  private:
    //! Don't permit copy constructor
    MKLRNG(const MKLRNG&) = delete;
    //! Don't permit copy assigment
    MKLRNG& operator=(const MKLRNG&) = delete;
    //! Don't permit move constructor
    MKLRNG(MKLRNG&&) = delete;
    //! Don't permit move assigment
    MKLRNG& operator=(MKLRNG&&) = delete;

    int m_nthreads;

    //! Random number stream for threads
    std::unique_ptr< VSLStreamStatePtr[] > m_stream;
};

} // quinoa::

#endif // MKLRNG_h
