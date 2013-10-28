//******************************************************************************
/*!
  \file      src/Paradigm/OpenMP.h
  \author    J. Bakosi
  \date      Sun 27 Oct 2013 03:32:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     OpenMP specifics
  \details   OpenMP specifics
*/
//******************************************************************************
#ifndef OpenMP_h
#define OpenMP_h

#include <Exception.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMP.h>

namespace tk {

//! OpenMP programming paradigm
class OpenMP {

  public:
    //! Constructor
    explicit OpenMP() :
      #ifdef _OPENMP
        m_available(true),
        m_used(true),
        m_nthreads(omp_get_max_threads())
      #else  // _OPENMP
        m_available(false),
        m_used(false),
        m_nthreads(1)
      #endif // _OPENMP
    { Assert(m_nthreads != 0, ExceptType::FATAL, "Need at least one thread"); }

    //! Destructor
    ~OpenMP() noexcept = default;

    //! Return true if compiled with OpenMP
    //! \return true if compiled with OpenMP enabled
    bool available() const noexcept { return m_available; }

    //! Query if OpenMP is used
    //! \return true if OpenMP is used
    bool used() const noexcept { return m_used; }

    //! Constant accessor to number of OpenMP threads
    int nthreads() const noexcept { return m_nthreads; }

  private:
    //! Don't permit copy constructor
    OpenMP(const OpenMP&) = delete;
    //! Don't permit copy assigment
    OpenMP& operator=(const OpenMP&) = delete;
    //! Don't permit move constructor
    OpenMP(OpenMP&&) = delete;
    //! Don't permit move assigment
    OpenMP& operator=(OpenMP&&) = delete;

    const bool m_available;           //!< True if OpenMP is available
    const bool m_used;                //!< True if OpenMP is used
    const int m_nthreads;             //!< Number of OpenMP threads
};

} // tk::

#endif // OpenMP_h
