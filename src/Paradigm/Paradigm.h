//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.h
  \author    J. Bakosi
  \date      Fri May 31 12:22:22 2013
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
    explicit Paradigm() = default;

    //! Destructor
    ~Paradigm() = default;

    //! Echo paradigm and configuration
    void echo() const;

    //! Query if OpenMP is available
    bool availOpenMP() const noexcept { return m_omp.available(); }

    //! Query if OpenMP is used
    bool usedOpenMP() const noexcept { return m_omp.used(); }

    //! Const accessor to OpenMP object
    const OpenMP* getOpenMP() const noexcept { return &m_omp; }

    //! Accessor to number of OpenMP threads
    int nthread() const noexcept { return m_omp.nthread(); }

  private:
    //! Don't permit copy constructor
    Paradigm(const Paradigm&) = delete;
    //! Don't permit copy assigment
    Paradigm& operator=(const Paradigm&) = delete;
    //! Don't permit move constructor
    Paradigm(Paradigm&&) = delete;
    //! Don't permit move assigment
    Paradigm& operator=(Paradigm&&) = delete;

    const OpenMP m_omp;       //!< OpenMP
};

} // namespace Quinoa

#endif // Paradigm_h
