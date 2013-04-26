//******************************************************************************
/*!
  \file      src/Random/MKLException.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:02:33 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKLException class declaration
  \details   MKLException class declaration
*/
//******************************************************************************
#ifndef MKLException_h
#define MKLException_h

#include <string>
#include <map>

using namespace std;

#include <RandomException.h>

namespace Quinoa {

//! MKL exception types
enum MKLExceptType { MKL_UNKNOWN_METHOD=0,
                     MKL_UNKNOWN_TABLE,
                     MKL_UNKNOWN_STREAM,
                     MKL_BAD_NTHREADS,
                     MKL_BAD_NUMBER,
                     MKL_VSL_ERROR,
                     NUM_MKL_EXCEPT
};

//! MKL exception error messages
const string MKLMsg[NUM_MKL_EXCEPT] = {
  "Unknown VSL generation method",
  "Random number table not found",
  "Random number stream not found",
  "Wrong number of threads",
  "Bad number of items"
  "VSL ",
};

//! MKLException : RandomException
class MKLException : public RandomException {

  public:
    //! Constructor
    explicit MKLException(const ExceptType except,
                          const MKLExceptType mklExcept,
                          const string& file,
                          const string& func,
                          const unsigned int& line) :
      RandomException(except, RND_MKL, file, func, line), m_except(mklExcept) {}

    //! Move constructor for throws, default compiler generated
    MKLException(MKLException&&) = default;

    //! Destructor
    virtual ~MKLException() noexcept = default;

    //! Handle MKLException
    virtual ErrCode handleException(Driver* const driver);

  protected:
    //! Permit copy constructor only for children
    MKLException(const MKLException&) = default;

  private:
    //! Don't permit copy assignment
    MKLException& operator=(const MKLException&) = delete;
    //! Don't permit move assignment
    MKLException& operator=(MKLException&&) = delete;

    //! MKL exception type (MKL_UNIMPLEMENTED, MKL_UNKNOWN_METHOD, etc.)
    const MKLExceptType m_except;
};

} // namespace Quinoa

#endif // MKLException_h
