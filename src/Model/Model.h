//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:20:47 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <string>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;

//! Model base
class Model {

  public:
    //! Constructor
    Model(Memory* const memory,
          Paradigm* const paradigm,
          Control* const control,
          const string& name);

    //! Destructor
    virtual ~Model() {};

    //! Name accessor
    const string& name() { return m_name; }

  protected:
    Memory* const m_memory;           //!< Memory object pointer
    Paradigm* const m_paradigm;       //!< Parallel programming object pointer
    Control* const m_control;         //!< Parallel programming object pointer

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;

    const string m_name;          //!< Name of model
};

} // namespace Quinoa

#endif // Model_h
