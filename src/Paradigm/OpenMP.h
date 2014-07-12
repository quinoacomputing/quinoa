//******************************************************************************
/*!
  \file      src/Paradigm/OpenMP.h
  \author    J. Bakosi
  \date      Sun 25 May 2014 06:11:02 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     OpenMP specifics
  \details   OpenMP specifics
*/
//******************************************************************************
#ifndef OpenMP_h
#define OpenMP_h

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tk {

//! OpenMP programming paradigm
class OpenMP {

  public:
    //! Constructor
    explicit OpenMP() = default;

    //! Destructor
    ~OpenMP() = default;

    //! Return true if compiled with OpenMP
    //! \return true if compiled with OpenMP enabled
    bool available() const noexcept {
      #ifdef _OPENMP
      return true;
      #else
      return false;
      #endif
    }

    //! Query if OpenMP is used
    //! \return true if OpenMP is used
    bool used() const noexcept {
      #ifdef _OPENMP
      return true;
      #else
      return false;
      #endif
    }

    //! Constant accessor to number of OpenMP threads
    //! \return the number of OpenMP threads used
    int nthreads() const noexcept {
      #ifdef _OPENMP
        return omp_get_max_threads();
      #else
        return 1;
      #endif
    }

  private:
    //! Don't permit copy constructor
    OpenMP(const OpenMP&) = delete;
    //! Don't permit copy assigment
    OpenMP& operator=(const OpenMP&) = delete;
    //! Don't permit move constructor
    OpenMP(OpenMP&&) = delete;
    //! Don't permit move assigment
    OpenMP& operator=(OpenMP&&) = delete;
};

} // tk::

#endif // OpenMP_h
