//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Fri May 10 17:22:21 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <QuinoaTypes.h>
#include <Exception.h>

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
                   const int npar) : m_memory(memory),
                                     m_paradigm(paradigm),
                                     m_control(control),
                                     m_npar(npar) {
      ErrChk(m_npar > 0, FATAL, "Wrong number of particles");
    }

    //! Destructor
    virtual ~Model() noexcept = default;

  protected:
    Memory* const m_memory;           //!< Memory object pointer
    Paradigm* const m_paradigm;       //!< Parallel programming object pointer
    Control* const m_control;         //!< Parallel programming object pointer
    const int m_npar;                 //!< Number of particles

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
