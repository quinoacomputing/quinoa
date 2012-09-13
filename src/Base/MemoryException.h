//******************************************************************************
/*!
  \file      src/Base/MemoryException.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 06:09:12 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryException class declaration
  \details   MemoryException class declaration
*/
//******************************************************************************
#ifndef MemoryException_h
#define MemoryException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! MemoryException types
enum MemoryExceptionType { BAD_ALLOC=0,  //!< std::bad_alloc
                           BAD_INSERT,   //!< unsuccessful STL::insert
                           BAD_NAME,     //!< non-unique variable name
                           EMPTY_STORE,  //!< memory store empty
                           NOT_FOUND,    //!< STL::find did not find entry
                           NOT_ERASED,   //!< STL::erase did not erase entry
                           UNDEFINED,    //!< memory entry not defined
                           NUM_MEMORY_EXCEPTIONS
};

//! MemoryException error messages
const string MemoryMessage[NUM_MEMORY_EXCEPTIONS] = {
  "Cannot allocate memory",
  "Cannot insert new memory entry",
  "Memory entry already exists",
  "Memory story is empty",
  "Memory entry not found",
  "Memory entry not erased",
  "Memory entry not defined"
};

//! MemoryException : Exception
class MemoryException : Exception {

  public:
    //! Constructor
    MemoryException(ExceptionType exception,
                    MemoryExceptionType memoryException) :
      Exception(exception), m_exception(memoryException) {}

    //! Copy constructor
    MemoryException(const MemoryException&);

    //! Destructor
    ~MemoryException() {}

    //! Handle MemoryException
    ErrorCode handleException(Driver* driver);

  private:
    //! Dont' permit assigment operator
    MemoryException& operator=(const MemoryException&);

    //! MemoryException (BAD_ALLOC, BAD_INSERT, BAD_NAME, etc.)
    MemoryExceptionType m_exception;
};

} // namespace Quinoa

#endif // MemoryException_h
