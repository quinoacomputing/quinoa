//******************************************************************************
/*!
  \file      src/Statistics/Distribution.h
  \author    J. Bakosi
  \date      Wed 24 Apr 2013 11:20:40 PM MDT
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

using namespace std;

namespace Quinoa {

//! Distribution estimator base
class Distribution {

  public:
    //! Constructor
    explicit Distribution() : m_nsample(0) {}

    //! Destructor
    virtual ~Distribution() {};

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
