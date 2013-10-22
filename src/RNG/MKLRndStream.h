//******************************************************************************
/*!
  \file      src/RNG/MKLRndStream.h
  \author    J. Bakosi
  \date      Tue Oct 22 15:43:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation from MKL streams
  \details   Streams are used to generate a few random numbers with no
             restrictions on the distribution parameters using leap-frogging
             between threads
*/
//******************************************************************************
#ifndef MKLRndStream_h
#define MKLRndStream_h

#include <MKL.h>

namespace tk {

//! MKL-based random number generator from leap-frog streams
class MKLRndStream : public MKL {

  public:
    //! Constructor: Create random number generator leap-frog stream
    explicit MKLRndStream(int nthread,
                          int brng,
                          unsigned int seed);

    //! Destructor: Destroy random number generator leap-frog stream
    ~MKLRndStream() noexcept override;

    //! Constant accessor to VSL streams
    const VSLStreamStatePtr* getStr() const noexcept { return m_stream; }

  private:
    //! Don't permit copy constructor
    MKLRndStream(const MKLRndStream&) = delete;
    //! Don't permit copy assigment
    MKLRndStream& operator=(const MKLRndStream&) = delete;
    //! Don't permit move constructor
    MKLRndStream(MKLRndStream&&) = delete;
    //! Don't permit move assigment
    MKLRndStream& operator=(MKLRndStream&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    const int m_nthread;             //!< Number of threads to use
    VSLStreamStatePtr* m_stream;     //!< Array of pointers to thread-streams
};

} // namespace tk

#endif // MKLRndStream_h
