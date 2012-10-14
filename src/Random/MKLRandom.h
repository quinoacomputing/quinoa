//******************************************************************************
/*!
  \file      src/Base/MKLRandom.h
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 10:43:41 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRandom_h
#define MKLRandom_h

#include <vector>

#include <mkl_vsl.h>

#include <Random.h>

using namespace std;

namespace Quinoa {

//! Probability distributions
enum Distribution { UNIFORM=0,        //!< Uniform
                    GAUSSIAN,         //!< Gaussian
                    GAMMA,            //!< Gamma
                    NUM_DIST_TYPES
};

//! MKL-based random number generator
class MKLRandom : Random {

  public:
    //! Add random table
    void addTable(Distribution dist, size_t number);

    //! Constructor: Setup random number generator streams
    MKLRandom(const Int nthreads, const uInt seed) : Random(nthreads, seed) {}

    //! Destructor: Destroy random number generator streams
    ~MKLRandom();

  private:
    //! Don't permit copy constructor
    MKLRandom(const MKLRandom&) = delete;
    //! Don't permit copy assigment
    MKLRandom& operator=(const MKLRandom&) = delete;
    //! Don't permit move constructor
    MKLRandom(MKLRandom&&) = delete;
    //! Don't permit move assigment
    MKLRandom& operator=(MKLRandom&&) = delete;

    //! Call MKL's vslNewStream() and handle error
    void newStream(VSLStreamStatePtr* stream,
                   const int& brng,
                   const unsigned int& seed);

    //! Stream tables to generate fixed numbers of random numbers with fixed
    //! properties using Random::m_nthreads
    vector<VSLStreamStatePtr*> table;
};

} // namespace Quinoa

#endif // MKLRandom_h
