//******************************************************************************
/*!
  \file      src/Paradigm/OpenMP.h
  \author    J. Bakosi
  \date      Fri May 31 12:21:59 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     OpenMP specifics
  \details   OpenMP specifics
*/
//******************************************************************************
#ifndef OpenMP_h
#define OpenMP_h

namespace Quinoa {

//! OpenMP programming paradigm
class OpenMP {

  public:
    //! Constructor
    explicit OpenMP();

    //! Destructor
    ~OpenMP() noexcept = default;

    //! Return true if compiled with OpenMP
    //! \return true if compiled with OpenMP enabled
    bool available() const noexcept { return m_available; }

    //! Query if OpenMP is used
    //! \return true if OpenMP is used
    bool used() const noexcept { return m_used; }

    //! Constant accessor to number of OpenMP threads
    int nthread() const noexcept { return m_nthread; }

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
    const int m_nthread;              //!< Number of OpenMP threads
};

} // namespace Quinoa

#endif // OpenMP_h
