//******************************************************************************
/*!
  \file      src/Base/Driver.h
  \author    J. Bakosi
  \date      Tue 04 Sep 2012 10:19:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class declaration
  \details   Driver base class declaration
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <Memory.h>

namespace Quinoa {

//! Driver base class declaration
class Driver {

  public:
    //! Constructor
    Driver(Memory* memory) : m_memory(memory) {}

    //! Destructor
    ~Driver();

    //! Finalize (either at normal exit, or due to exception)
    void finalize();

  private:
    //! Don't permit copy operator
    Driver(const Driver&);
    //! Don't permit assigment operator
    Driver& operator=(const Driver&);

    //! Pointer to Memory object
    Memory* m_memory;
};

} // namespace Quinoa

#endif // Driver_h
