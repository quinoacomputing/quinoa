//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.h
  \author    J. Bakosi
  \date      Tue 06 Nov 2012 06:08:06 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************
#ifndef Paradigm_h
#define Paradigm_h

#include <OpenMP.h>

namespace Quinoa {

//! Parallel programming paradigms
class Paradigm {

  public:
    //! Constructor
    Paradigm() = default;

    //! Destructor
    ~Paradigm() = default;

    //! Echo paradigm and configuration
    void echo();

    //! Query if OpenMP is available
    bool availOpenMP() const { return m_omp.available(); }

    //! Query if OpenMP is used
    bool usedOpenMP() const { return m_omp.used(); }

    //! Const accessor to OpenMP object
    const OpenMP* getOpenMP() const { return &m_omp; }

  private:
    //! Don't permit copy constructor
    Paradigm(const Paradigm&) = delete;
    //! Don't permit copy assigment
    Paradigm& operator=(const Paradigm&) = delete;
    //! Don't permit move constructor
    Paradigm(Paradigm&&) = delete;
    //! Don't permit move assigment
    Paradigm& operator=(Paradigm&&) = delete;

    OpenMP m_omp;       //!< OpenMP
};

} // namespace Quinoa

#endif // Paradigm_h
