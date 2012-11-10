//******************************************************************************
/*!
  \file      src/Base/MemoryException.h
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 09:25:49 AM MST
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

//! Memory exception types
enum MemExceptType { BAD_ALLOC=0,     //!< std::bad_alloc
                     BAD_INSERT,      //!< unsuccessful STL::insert
                     NONUNIQUE_NAME,  //!< non-unique variable name
                     EMPTY_STORE,     //!< memory store empty
                     NOT_FOUND,       //!< STL::find did not find entry
                     NOT_ERASED,      //!< STL::erase did not erase entry
                     UNDEFINED,       //!< memory entry not defined
                     BAD_NUMBER,      //!< bad number of items
                     EMPTY_NAME,      //!< bad number of items
                     NUM_MEM_EXCEPT
};

//! Memory exception error messages
const string MemMsg[NUM_MEM_EXCEPT] = {
  "Cannot allocate memory",
  "Cannot insert new memory entry",
  "Memory entry already exists",
  "Memory story is empty",
  "Memory entry not found",
  "Memory entry not erased",
  "Memory entry not defined",
  "Bad number of items"
  "Must specify a name with non-zero length"
};

//! MemoryException : Exception
class MemoryException : public Exception {

  public:
    //! Constructor
    MemoryException(ExceptType except, MemExceptType memExcept) :
      Exception(except), m_except(memExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    MemoryException(MemoryException&&) = default;

    //! Destructor
    virtual ~MemoryException() {}

    //! Handle MemoryException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy constructor
    MemoryException(const MemoryException&) = delete;
    //! Don't permit copy assignment
    MemoryException& operator=(const MemoryException&) = delete;
    //! Don't permit move assignment
    MemoryException& operator=(MemoryException&&) = delete;

    //! Memory exception type (BAD_ALLOC, BAD_INSERT, BAD_NAME, etc.)
    MemExceptType m_except;
};

} // namespace Quinoa

#endif // MemoryException_h
