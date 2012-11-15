//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Thu Nov 15 15:02:30 2012
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

//! Model base
class Model {

  public:
    //! Constructor
    Model(Memory* memory, Paradigm* paradigm, const string& name);

    //! Destructor
    virtual ~Model();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

  protected:
    Memory* m_memory;             //!< Memory object pointer
    Paradigm* m_paradigm;         //!< Parallel programming object pointer
    const string m_name;          //!< Name of model

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;
};

} // namespace Quinoa

#endif // Model_h
