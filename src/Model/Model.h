//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Fri 16 Nov 2012 09:48:04 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <string>
#include <limits>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;

//! Model base
class Model {

  public:
    //! Constructor
    Model(Memory* memory,
          Paradigm* paradigm,
          const string& name,
          const real time,
          const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~Model();

    //! Echo informaion on model
    virtual void echo() = 0;

    //! Initialize model
    virtual void init() = 0;

    //! Solve model
    virtual void solve() = 0;

  protected:
    Memory* m_memory;             //!< Memory object pointer
    Paradigm* m_paradigm;         //!< Parallel programming object pointer
    const string m_name;          //!< Name of model
    const real m_time;            //!< Maximum time to simulate
    const int m_nstep;            //!< Maximum number of time steps to take

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
