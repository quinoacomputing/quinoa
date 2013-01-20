//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 01:03:44 PM MST
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

    //! One-liner report
    virtual void report(const int it,
                        const real t,
                        const real dt,
                        long int& hrs2beg,
                        long int& mins2beg,
                        long int& secs2beg,
                        long int& hrs2end,
                        long int& mins2end,
                        long int& secs2end);

  protected:
    Memory* m_memory;             //!< Memory object pointer
    Paradigm* m_paradigm;         //!< Parallel programming object pointer
    const string m_name;          //!< Name of physics
    const real m_time;            //!< Maximum time to simulate
    const int m_echo;             //!< One-liner info in every few time step
    const int m_nstep;            //!< Maximum number of time steps to take
    struct timeval m_startTime;   //!< Date/time when time-adv loop started

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
