//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:56:46 AM KST
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
enum ValueType { BOOL=0,        //!< Boolean value
                 INT,           //!< Integer value
                 REAL,          //!< Real value
                 NUM_VALUE_TYPES
};

//! Size of value types
const size_t SizeOf[NUM_VALUE_TYPES] = { sizeof(Bool),  //!< Size of Bool
                                         sizeof(Int),   //!< Size of Integer
                                         sizeof(Real)   //!< Size of Real
};

//! Name of value types
const string ValueName[NUM_VALUE_TYPES] = { "bool",  //! Screen name of bool
                                            "int",   //! Screen name of integer
                                            "real"   //! Screen name of real
};

//! Variable types
enum VariableType { SCALAR=0,       //!< Scalar quantity
                    VECTOR,         //!< Vector quantity
                    SYMTENSOR,      //!< Symmetric tensor quantity
                    TENSOR,         //!< Tensor quantity
                    NUM_VARIABLE_TYPES
};

//! Variable components
const Int VariableComponents[NUM_VARIABLE_TYPES] { 1,  //!< Scalar
                                                   3,  //!< Vector
                                                   6,  //!< Symmetric tensor
                                                   9   //!< Tensor
};

//! Name of variable types
const string VariableTypeName[NUM_VARIABLE_TYPES] = { "scalar",
                                                      "vector",
                                                      "symtensor",
                                                      "tensor"
};

//! Output width of MemoryEntry fields
const Int EntryWidth[] = { 10,  //! Width of Name field
                           10,  //! Width of Number field
                            5,  //! Width of Value field
                            9,  //! Width of Value size field
                           10,  //! Width of Variable field
                           10,  //! Width of Bytes field
                            6,  //! Width of Plot field
                            7,  //! Width of Restart field
                           10   //! Width of Ptr field
};

//! MemoryEntry base class
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
      if (m_ptr) delete [] static_cast<char*>(m_ptr);
    }

    //! Don't permit copy operator
    MemoryEntry(const MemoryEntry&);
    //! Don't permit assigment operator
    MemoryEntry& operator=(const MemoryEntry&);

    //! One-liner accessor for all fields
    string line();

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
