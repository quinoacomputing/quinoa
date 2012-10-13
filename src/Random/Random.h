//******************************************************************************
/*!
  \file      src/Base/Random.h
  \author    J. Bakosi
  \date      Thu 11 Oct 2012 08:29:27 PM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef Random_h
#define Random_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! Random number generator base
class Random {

  protected:
    //! Constructor: Setup random number generators
    Random(const Int nthreads, const uInt seed) :
      m_nthreads(nthreads), m_seed(seed) {}

    //! Destructor: Destroy random number generators
    ~Random() = default;

    Int m_nthreads;              //!< Number of threads
    const uInt m_seed;           //!< Seed

  private:
    //! Don't permit copy constructor
    Random(const Random&) = delete;
    //! Don't permit copy assigment
    Random& operator=(const Random&) = delete;
    //! Don't permit move constructor
    Random(Random&&) = delete;
    //! Don't permit move assigment
    Random& operator=(Random&&) = delete;
};

} // namespace Quinoa

#endif // Random_h
