//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Fri Apr 26 16:21:30 2013
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
    explicit Model(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   const string& name);

    //! Destructor
    virtual ~Model() noexcept = default;

    //! Name accessor
    const string& name() const { return m_name; }

    //! Interface to accessor to number of particle properties
    virtual int nprop() const = 0;

    //! Constant accessor to particle properties pointer
    virtual const real* particles() const = 0;

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
