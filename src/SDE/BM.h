//******************************************************************************
/*!
  \file      src/SDE/BM.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:38:45 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************
#ifndef BM_h
#define BM_h

#include <SDE.h>

namespace quinoa {

//! BM : Mass
class BM : public Mass {

  public:
    //! Constructor
    explicit BM() = default;

    //! Destructor
    ~BM() noexcept override = default;

  private:
    //! Don't permit copy constructor
    BM(const BM&) = delete;
    //! Don't permit copy assigment
    BM& operator=(const BM&) = delete;
    //! Don't permit move constructor
    BM(BM&&) = delete;
    //! Don't permit move assigment
    BM& operator=(BM&&) = delete;
};

} // quinoa::

#endif // BM_h
