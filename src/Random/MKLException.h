//******************************************************************************
/*!
  \file      src/Random/MKLException.h
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 05:06:15 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKLException class declaration
  \details   MKLException class declaration
*/
//******************************************************************************
#ifndef MKLException_h
#define MKLException_h

#include <string>

using namespace std;

#include <RandomException.h>

namespace Quinoa {

//! MKL exception types
enum MKLExceptType { MKL_UNIMPLEMENTED=0,  //!< VSL feature not implemented
                     NUM_MKL_EXCEPT
};

//! MKL exception error messages
const string MKLMsg[NUM_MKL_EXCEPT] = {
  "VSL feature not yet implemented"
};

//! MKLException : RandomException
class MKLException : protected RandomException {

  public:
    //! Constructor
    MKLException(ExceptType except, MKLExceptType mklExcept) :
      RandomException(except, RND_MKL), m_except(mklExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    MKLException(MKLException&&) = default;

    //! Destructor
    ~MKLException() = default;

    //! Handle MKLException
    ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy constructor
    MKLException(const MKLException&) = delete;
    //! Don't permit copy assignment
    MKLException& operator=(const MKLException&) = delete;
    //! Don't permit move assignment
    MKLException& operator=(MKLException&&) = delete;

    //! MKL exception type (, etc.)
    MKLExceptType m_except;
};

} // namespace Quinoa

#endif // MKLException_h
