//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Tue May  7 10:44:49 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;
class Timer;

//! Physics base
class Physics {

  public:
    //! Constructor
    explicit Physics(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     Timer* const timer) noexcept;

    //! Destructor
    virtual ~Physics() noexcept = default;

    //! Echo informaion on physics
    virtual void echo() const = 0;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

  protected:
    Memory* const m_memory;       //!< Memory object
    Paradigm* const m_paradigm;   //!< Parallel programming object
    Control* const m_control;     //!< Control object
    Timer* const m_timer;         //!< Timer object
    const int m_nthread;          //!< Number of threads

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
