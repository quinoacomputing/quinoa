//******************************************************************************
/*!
  \file      src/Models/Mix/MixModel.h
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 08:23:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MixModel base
  \details   MixModel base
*/
//******************************************************************************
#ifndef MixModel_h
#define MixModel_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! MixModel base
class MixModel {

  public:
    //! Constructor
    MixModel() = default;

    //! Destructor
    ~MixModel() = default;

    //! Interface for setting initial conditions
    virtual void setIC() = 0;

  private:
    //! Don't permit copy constructor
    MixModel(const MixModel&) = delete;
    //! Don't permit copy assigment
    MixModel& operator=(const MixModel&) = delete;
    //! Don't permit move constructor
    MixModel(MixModel&&) = delete;
    //! Don't permit move assigment
    MixModel& operator=(MixModel&&) = delete;
};

} // namespace Quinoa

#endif // MixModel_h
