//******************************************************************************
/*!
  \file      src/Geometry/Geometry.h
  \author    J. Bakosi
  \date      Thu Jun 20 06:48:26 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************
#ifndef Geometry_h
#define Geometry_h

#include <QuinoaTypes.h>
#include <Exception.h>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;

//! Geometry base
class Geometry {

  public:
    //! Constructor
    explicit Geometry(Memory* const memory,
                      Paradigm* const paradigm,
                      Control* const control,
                      Timer* const timer) :
      m_memory(memory),
      m_paradigm(paradigm),
      m_control(control),
      m_timer(timer) {}

    //! Destructor
    virtual ~Geometry() noexcept {}

    //! Initialize geometry
    virtual void init() = 0;

    //! Space-fill geometry
    virtual void fill() = 0;

  protected:
    Memory* const m_memory;                    //!< Memory object
    Paradigm* const m_paradigm;                //!< Parallel programming object
    Control* const m_control;                  //!< Control object
    Timer* const m_timer;                      //!< Timer object

  private:
    //! Don't permit copy constructor
    Geometry(const Geometry&) = delete;
    //! Don't permit copy assigment
    Geometry& operator=(const Geometry&) = delete;
    //! Don't permit move constructor
    Geometry(Geometry&&) = delete;
    //! Don't permit move assigment
    Geometry& operator=(Geometry&&) = delete;
};

} // namespace Quinoa

#endif // Geometry_h
