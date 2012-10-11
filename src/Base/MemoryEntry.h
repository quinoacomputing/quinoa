//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.h
  \author    J. Bakosi
  \date      Wed 10 Oct 2012 01:46:10 PM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory entry
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
enum ValType { BOOL=0,        //!< Boolean value
               INT,           //!< Integer value
               REAL,          //!< Real value
               NUM_VAL_TYPES
};
//! Size of value types
constexpr size_t SizeOf[NUM_VAL_TYPES] = { sizeof(Bool),
                                           sizeof(Int),
                                           sizeof(Real)
};
//! (Screen) names of value types
const string ValName[NUM_VAL_TYPES] = { "Bool",
                                        "Int",
                                        "Real"
};

//! Variable types
enum VarType { SCALAR=0,   //!< Scalar quantity
               VECTOR,     //!< Vector quantity
               SYMTENSOR,  //!< Symmetric tensor quantity
               TENSOR,     //!< Tensor quantity
               NUM_VAR_TYPES
};
//! Variable components
const Int VarComp[NUM_VAR_TYPES] { 1,  //!< Scalar
                                   3,  //!< Vector
                                   6,  //!< Symmetric tensor
                                   9   //!< Tensor
};
//! Name of variable types
const string VarTypeName[NUM_VAR_TYPES] = { "scalar",
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

//! Memory entry field designators
enum MemoryEntryField { UNSPECIFIED=0,
                        BYTES,
                        NUMBER,
                        VALUE,
                        VARIABLE,
                        NAME,
                        PLOT,
                        RESTART,
                        PTR
};

//! Memory entry
class MemoryEntry {

  private:
    //! Befriend class Memory to allow direct manipulation of private fields
    friend class Memory;

    //! Constructor: fill in all fields
    MemoryEntry(size_t bytes,
                size_t number,
                ValType value,
                VarType variable,
                string name,
                Bool plot,
                Bool restart,
                void* ptr) :
      m_bytes(bytes),
      m_number(number),
      m_value(value),
      m_variable(variable),
      m_name(name),
      m_plot(plot), 
      m_restart(restart),
      m_ptr(ptr) {}

    //! Destructor: free allocated memory when leaving scope
    ~MemoryEntry() {
      if (m_ptr) delete [] static_cast<char*>(m_ptr);
    }

    //! Don't permit copy constructor
    MemoryEntry(const MemoryEntry&) = delete;
    //! Don't permit copy assigment
    MemoryEntry& operator=(const MemoryEntry&) = delete;
    //! Don't permit move constructor
    MemoryEntry(MemoryEntry&&) = delete;
    //! Don't permit move assigment
    MemoryEntry& operator=(MemoryEntry&&) = delete;

    //! One-liner accessor for all fields
    string line();
 
    size_t m_bytes;           //!< Size in bytes (number of chars) allocated
    size_t m_number;          //!< Number of items
    ValType m_value;          //!< Value type (BOOL, INT, etc.)
    VarType m_variable;       //!< Variable type (SCALAR, VECTOR, etc.)
    string m_name;            //!< Variable name
    Bool m_plot;              //!< Variable can be plotted
    Bool m_restart;           //!< Write to restart file
    void* m_ptr;              //!< Pointer to data
};

} // namespace Quinoa

#endif // MemoryEntry_h
