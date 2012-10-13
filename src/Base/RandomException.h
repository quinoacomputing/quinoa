//******************************************************************************
/*!
  \file      src/Base/RandomException.h
  \author    J. Bakosi
  \date      Thu 11 Oct 2012 08:47:38 PM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RandomException class declaration
  \details   RandomException class declaration
*/
//******************************************************************************
#ifndef RandomException_h
#define RandomException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! RandomException types
enum RndExceptType { RND_UNIMPLEMENTED=0,   //!< unimplemented feature
                     NUM_RND_EXCEPT
};

//! Random exception error messages
const string RndMsg[NUM_RND_EXCEPT] = {
  "MKL VSL feature not yet implemented",
};

//! RandomException : Exception
class RandomException : Exception {

  public:
    //! Constructor with message from thrower
    RandomException(ExceptType except, RndExceptType rndExcept) :
      Exception(except), m_except(rndExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    RandomException(RandomException&&) = default;

    //! Destructor
    ~RandomException() = default;

    //! Handle RandomException
    ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy constructor
    RandomException(const RandomException&) = delete;
    //! Don't permit copy assignment
    RandomException& operator=(const RandomException&) = delete;
    //! Don't permit move assignment
    RandomException& operator=(RandomException&&) = delete;

    //! Random exception type (RND_UNIMPLEMENTED, etc.)
    RndExceptType m_except;
};

} // namespace Quinoa

#endif // RandomException_h
