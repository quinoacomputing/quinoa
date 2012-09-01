//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.h
  \author    J. Bakosi
  \date      Sat 01 Sep 2012 03:04:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory entry declaration
  \details   The memory store contains memory entries
*/
//******************************************************************************
#ifndef MemoryEntry_h
#define MemoryEntry_h

#include <string>
#include <map>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! Variable types
enum VarType { SCALAR_VAR = 0,     //!< Scalar quantity
               VECTOR_VAR,         //!< Vector quantity
               SYMTENSOR_VAR,      //!< Symmetric tensor quantity
               TENSOR_VAR,         //!< Tensor quantity
               NUM_VAR_TYPES
};

//! Variable components
const Int VarComp[NUM_VAR_TYPES] { 1,  //!< Scalar
                                   3,  //!< Vector
                                   6,  //!< Symmetric tensor
                                   9   //!< Tensor
};

//! MemoryEntry declaration
template<class V>
class MemoryEntry {

  //! Befriend class Memory to allow direct manipulation of private fields
  template<class> friend class Memory;

  private:
    //! Constructor
    MemoryEntry(size_t number,
                VarType type,
                string name,
                Bool plot,
                Bool restart,
                V* ptr) :
      m_number(number),
      m_type(type),
      m_name(name),
      m_plot(plot), 
      m_restart(restart),
      m_ptr(ptr) {}

    //! Destructor: Free allocated memory when leaving scope
    ~MemoryEntry() {
      delete [] m_ptr;
    }

    //! Don't permit copy operator
    MemoryEntry(const MemoryEntry&);
    //! Don't permit assigment operator
    MemoryEntry& operator=(const MemoryEntry&);

    size_t m_number;        //!< Number of items
    VarType m_type;         //!< Variable type (SCALAR_VAR, VECTOR_VAR, etc.)
    string m_name;          //!< Variable name
    Bool m_plot;            //!< Variable can be plotted
    Bool m_restart;         //!< Write to restart file
    V* m_ptr;               //!< Pointer to data
};

} // namespace Quinoa

#endif // MemoryEntry_h
