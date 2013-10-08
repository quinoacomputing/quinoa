//******************************************************************************
/*!
  \file      src/Statistics/Distribution.h
  \author    J. Bakosi
  \date      Mon 07 Oct 2013 08:44:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Distribution estimator base
  \details   Distribution estimator base
*/
//******************************************************************************
#ifndef Distribution_h
#define Distribution_h

namespace tk {

//! Distribution estimator base
class Distribution {

  public:
    //! Constant accessor to PDF map
    virtual const int& getNsample() const = 0;

  protected:
    //! Constructor: designed to be base-only
    explicit Distribution() noexcept : m_nsample(0) {}

    //! Destructor: designed to be freed via children-only
    virtual ~Distribution() noexcept = default;

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

} // namespace tk

#endif // Distribution_h
