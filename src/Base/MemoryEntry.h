//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.h
  \author    J. Bakosi
  \date      Wed 05 Sep 2012 08:28:58 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory entry declaration
  \details   The memory store contains memory entries
*/
//******************************************************************************
#ifndef MemoryEntry_h
#define MemoryEntry_h

#include <string>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! Value types
enum ValueType { BOOL_VAL=0,        //!< Boolean value
                 INT_VAL,           //!< Integer value
                 REAL_VAL,          //!< Real value
                 NUM_VALUE_TYPES
};

//! Size of value types
const size_t SizeOf[NUM_VALUE_TYPES] = { sizeof(Bool),  //!< Size of Bool
                                         sizeof(Int),   //!< Size of Integer
                                         sizeof(Real)   //!< Size of Real
};


//! Variable types
enum VariableType { SCALAR_VAR=0,       //!< Scalar quantity
                    VECTOR_VAR,         //!< Vector quantity
                    SYMTENSOR_VAR,      //!< Symmetric tensor quantity
                    TENSOR_VAR,         //!< Tensor quantity
                    NUM_VARIABLE_TYPES
};

//! Variable components
const Int VariableComponents[NUM_VARIABLE_TYPES] { 1,  //!< Scalar
                                                   3,  //!< Vector
                                                   6,  //!< Symmetric tensor
                                                   9   //!< Tensor
};

//! MemoryEntry declaration
class MemoryEntry {

  //! Befriend class Memory to allow direct manipulation of private fields
  friend class Memory;

  private:
    //! Constructor
    MemoryEntry(size_t nbytes,
                size_t number,
                ValueType value,
                VariableType variable,
                string name,
                Bool plot,
                Bool restart,
                void* ptr) :
      m_nbytes(nbytes),
      m_number(number),
      m_value(value),
      m_variable(variable),
      m_name(name),
      m_plot(plot), 
      m_restart(restart),
      m_ptr(ptr) {}

    //! Destructor: Free allocated memory when leaving scope
    ~MemoryEntry() {
      delete [] static_cast<char*>(m_ptr);
    }

    //! Don't permit copy operator
    MemoryEntry(const MemoryEntry&);
    //! Don't permit assigment operator
    MemoryEntry& operator=(const MemoryEntry&);

    size_t m_nbytes;          //!< Size in bytes (number of chars) allocated
    size_t m_number;          //!< Number of items
    ValueType m_value;        //!< Value type (BOOL_VAL, INT_VAL, etc.)
    VariableType m_variable;  //!< Variable type (SCALAR_VAR, VECTOR_VAR, etc.)
    string m_name;            //!< Variable name
    Bool m_plot;              //!< Variable can be plotted
    Bool m_restart;           //!< Write to restart file
    void* m_ptr;              //!< Pointer to data
};

} // namespace Quinoa

#endif // MemoryEntry_h
