//******************************************************************************
/*!
  \file      src/Base/MemoryException.h
  \author    J. Bakosi
  \date      Fri Sep 14 17:25:47 2012
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
enum class MemExceptType : Int {
                           BAD_ALLOC=0,  //!< std::bad_alloc
                           BAD_INSERT,   //!< unsuccessful STL::insert
                           BAD_NAME,     //!< non-unique variable name
                           EMPTY_STORE,  //!< memory store empty
                           NOT_FOUND,    //!< STL::find did not find entry
                           NOT_ERASED,   //!< STL::erase did not erase entry
                           UNDEFINED,    //!< memory entry not defined
                           NUM_MEM_EXCEPT
};
//! Number of memory exception types
const Int NUM_MEM_EXCEPT = static_cast<Int>(MemExceptType::NUM_MEM_EXCEPT);

//! Memory exception error messages
const string MemMsg[NUM_MEM_EXCEPT] = {
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
    MemoryException(ExceptType except, MemExceptType memExcept) :
      Exception(except), m_except(memExcept) {}

    //! Copy constructor
    MemoryException(const MemoryException&);

    //! Destructor
    ~MemoryException() {}

    //! Handle MemoryException
    ErrCode handleException(Driver* driver);

  private:
    //! Dont' permit assigment operator
    MemoryException& operator=(const MemoryException&);

    //! Memory exception type (BAD_ALLOC, BAD_INSERT, BAD_NAME, etc.)
    MemExceptType m_except;
};

} // namespace Quinoa

#endif // MemoryException_h
