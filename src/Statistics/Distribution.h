//******************************************************************************
/*!
  \file      src/Statistics/Distribution.h
  \author    J. Bakosi
  \date      Tue Jul  2 16:11:20 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Distribution estimator base
  \details   Distribution estimator base
*/
//******************************************************************************
#ifndef Distribution_h
#define Distribution_h

#include <vector>
#include <iostream>

#include <QuinoaTypes.h>

namespace Quinoa {

//! Distribution estimator base
class Distribution {

  public:
    //! Constructor
    explicit Distribution() noexcept : m_nsample(0) {}

    //! Destructor
    virtual ~Distribution() noexcept = default;

    //! Constant accessor to PDF map
    virtual const int& getNsample() const = 0;

  protected:
    int m_nsample;          //!< Number of samples collected

  private:
    //! Don't permit copy constructor
    Distribution(const Distribution&) = delete;
    //! Don't permit copy assigment
    Distribution& operator=(const Distribution&) = delete;
    //! Don't permit move constructor
    Distribution(Distribution&&) = delete;
    //! Don't permit move assigment
    Distribution& operator=(Distribution&&) = delete;
};

} // namespace Quinoa

#endif // Distribution_h
