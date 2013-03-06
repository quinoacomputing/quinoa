//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.h
  \author    J. Bakosi
  \date      Wed 06 Mar 2013 06:38:44 AM MST
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

    //! Accessor to number of OpenMP threads
    int nthread() const { return m_omp.nthread(); }

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
