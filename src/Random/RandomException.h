//******************************************************************************
/*!
  \file      src/Random/RandomException.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 06:19:56 PM MST
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
enum RndExceptType { RND_MKL=0,             //!< MKL error
                     RND_UNIMPLEMENTED,     //!< Unimplemented exception
                     NUM_RND_EXCEPT
};

//! Random exception error messages
const string RndMsg[NUM_RND_EXCEPT] = {
  "MKL exception: ",
  "Unimplemented random number generator exception"
};

//! RandomException : Exception
class RandomException : public Exception {

  public:
    //! Constructor with message from thrower
    RandomException(ExceptType except, RndExceptType rndExcept) :
      Exception(except), m_except(rndExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    //! Can only be thrown from within derived RandomException classes
    RandomException(RandomException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    RandomException(const RandomException&);

    //! Destructor
    virtual ~RandomException() {}

    //! Handle RandomException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    RandomException& operator=(const RandomException&) = delete;
    //! Don't permit move assignment
    RandomException& operator=(RandomException&&) = delete;

    //! Random exception type (RND_UNIMPLEMENTED, etc.)
    RndExceptType m_except;
};

} // namespace Quinoa

#endif // RandomException_h
