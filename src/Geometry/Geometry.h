//******************************************************************************
/*!
  \file      src/Geometry/Geometry.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:05:00 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************
#ifndef Geometry_h
#define Geometry_h

namespace quinoa {

class Memory;
class Paradigm;
class Timer;
class QuinoaControl;

//! Geometry base
class Geometry {

  public:
    //! Constructor
    explicit Geometry(Memory* const memory,
                      Paradigm* const paradigm,
                      QuinoaControl* const control,
                      Timer* const timer) noexcept :
      m_memory(memory),
      m_paradigm(paradigm),
      m_control(control),
      m_timer(timer) {}

    //! Destructor
    virtual ~Geometry() noexcept = default;

    //! Initialize geometry
    virtual void init() = 0;

    //! Space-fill geometry
    virtual void fill() = 0;

  protected:
    Memory* const m_memory;                    //!< Memory object
    Paradigm* const m_paradigm;                //!< Parallel programming object
    QuinoaControl* const m_control;            //!< Control object
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

} // namespace quinoa

#endif // Geometry_h
