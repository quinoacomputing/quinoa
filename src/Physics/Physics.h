//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Sat 19 Jan 2013 11:24:58 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <string>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;

//! Physics base
class Physics {

  public:
    //! Constructor
    Physics(Memory* memory,
            Paradigm* paradigm,
            const string& name,
            const real time,
            const int echo,
            const int nstep);

    //! Destructor
    virtual ~Physics();

    //! Echo informaion on physics
    virtual void echo() = 0;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

  protected:
    Memory* m_memory;             //!< Memory object pointer
    Paradigm* m_paradigm;         //!< Parallel programming object pointer
    const string m_name;          //!< Name of physics
    const real m_time;            //!< Maximum time to simulate
    const int m_echo;             //!< One-line info in every few time step
    const int m_nstep;            //!< Maximum number of time steps to take

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;
};

} // namespace Quinoa

#endif // Physics_h
