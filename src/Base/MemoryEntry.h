//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.h
  \author    J. Bakosi
  \date      Thu 30 Aug 2012 11:26:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory entry declaration
  \details   The memory store contains memory entries
*/
//******************************************************************************
#ifndef MemoryEntry_h
#define MemoryEntry_h

#include <string>
#include <map>

namespace Quinoa {

//! Variable types
enum VariableType { SCALAR_VAR = 0,     //!< Scalar quantity
                    VECTOR_VAR,         //!< Vector quantity
                    SYMTENSOR_VAR,      //!< Symmetric tensor quantity
                    TENSOR_VAR          //!< Tensor quantity
};
static pair<VariableType, size_t> VariableTypePair[] = {
  pair<VariableType, size_t>(SCALAR_VAR,      sizeof(double)),
  pair<VariableType, size_t>(VECTOR_VAR,    3*sizeof(double)),
  pair<VariableType, size_t>(SYMTENSOR_VAR, 6*sizeof(double)),
  pair<VariableType, size_t>(TENSOR_VAR,    9*sizeof(double)),
};
// Initialize a static const map based on the VariableTypePair array
static const map<VariableType, size_t>
  VariableSize(VariableTypePair,
               VariableTypePair +
               sizeof(VariableTypePair)/sizeof(VariableTypePair[0]));

//! Memory storage layout for multi-dimensional arrays
enum StorageLayout { LINEAR_ARRAY = 0,   //!< Linear memory storage
                     COLUMN_MAJOR_ARRAY, //!< Column-major storage
                     ROW_MAJOR_ARRAY,    //!< Row-major storage
                     NUM_MEMORY_LAYOUT
};

//! MemoryEntry declaration
class MemoryEntry {

  public:
    //! Constructor
    MemoryEntry(size_t number,
                VariableType type,
                int dim,
                string name,
                StorageLayout layout,
                bool plot,
                bool restart,
                void* ptr) :
      m_number(number),
      m_type(type),
      m_dim(dim),
      m_name(name),
      m_layout(layout),
      m_plot(plot), 
      m_restart(restart),
      m_ptr(ptr) {}

   //! Destructor
   ~MemoryEntry() {}

  private:
    //! Don't permit copy operator
    MemoryEntry(const MemoryEntry&);
    //! Don't permit assigment operator
    MemoryEntry& operator=(const MemoryEntry&);

    size_t m_number;            //!< Number of items per dimension
    VariableType m_type;        //!< Variable type
    int m_dim;                  //!< Number of dimensions
    string m_name;              //!< Variable name
    StorageLayout m_layout;     //!< Memory storage layout
    bool m_plot;                //!< Variable can be plotted
    bool m_restart;             //!< Write to restart file
    void* m_ptr;                //!< Pointer to data
};

} // namespace Quinoa

#endif // MemoryEntry.h
