//******************************************************************************
/*!
  \file      src/Random/MKLException.h
  \author    J. Bakosi
  \date      Mon Apr 29 15:27:51 2013
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
  "Bad number of items",
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
                          unsigned int line,
                          const string& message = "") noexcept :
      RandomException(except,
                      RND_MKL,
                      file,
                      func,
                      line,
                      MKLMsg[static_cast<int>(mklExcept)] + message) {}

    //! Move constructor for throws, default compiler generated
    MKLException(MKLException&&) = default;

    //! Destructor
    virtual ~MKLException() noexcept = default;

  protected:
    //! Permit copy constructor only for children
    MKLException(const MKLException&) = default;

  private:
    //! Don't permit copy assignment
    MKLException& operator=(const MKLException&) = delete;
    //! Don't permit move assignment
    MKLException& operator=(MKLException&&) = delete;
};

} // namespace Quinoa

#endif // MKLException_h
