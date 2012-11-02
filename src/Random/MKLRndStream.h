//******************************************************************************
/*!
  \file      src/Random/MKLRndStream.h
  \author    J. Bakosi
  \date      Thu 01 Nov 2012 07:31:17 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation from MKL streams
  \details   Streams are used to generate a few random numbers with no
             restrictions on the distribution parameters using leap-frogging
             between threads
*/
//******************************************************************************
#ifndef MKLRndStream_h
#define MKLRndStream_h

#include <QuinoaTypes.h>
#include <Memory.h>
#include <MKL.h>

using namespace std;

namespace Quinoa {

//! MKL-based random number generator from leap-frog streams
class MKLRndStream : public MKL {

  public:
    //! Constructor: Create random number generator leap-frog stream
    MKLRndStream(const int nthreads,
                 const int brng,
                 const unsigned int seed);

    //! Destructor: Destroy random number generator leap-frog stream
    ~MKLRndStream();

    //! Constant accessor to VSL streams
    const VSLStreamStatePtr* getStr() { return m_stream; }

  private:
    //! Don't permit copy constructor
    MKLRndStream(const MKLRndStream&) = delete;
    //! Don't permit copy assigment
    MKLRndStream& operator=(const MKLRndStream&) = delete;
    //! Don't permit move constructor
    MKLRndStream(MKLRndStream&&) = delete;
    //! Don't permit move assigment
    MKLRndStream& operator=(MKLRndStream&&) = delete;

    const int m_nthread;             //!< Number of threads to use
    VSLStreamStatePtr* m_stream;     //! Array of pointers to thread-streams
};

} // namespace Quinoa

#endif // MKLRndStream_h
