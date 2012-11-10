//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 06:56:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class declaration
  \details   Driver base class declaration
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <Memory.h>
#include <Model.h>

namespace Quinoa {

//! Driver base class
class Driver {

  public:
    //! Constructor
    Driver(Memory* memory);

    //! Destructor
    ~Driver();

    //! Setup
    void setup();

    //! Solve
    void solve();

    //! Finalize (either at normal exit, or due to exception)
    void finalize();

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;

    Memory* m_memory;           //!< Pointer to Memory object
    Model* m_model;             //!< Pointer to Model object
};

} // namespace Quinoa

#endif // Driver_h
