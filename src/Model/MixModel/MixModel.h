//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModel.h
  \author    J. Bakosi
  \date      Fri Nov 16 07:33:01 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MixModel base
  \details   MixModel base
*/
//******************************************************************************
#ifndef MixModel_h
#define MixModel_h

#include <string>

#include <QuinoaTypes.h>

namespace Quinoa {

using namespace std;

//! MixModel base
class MixModel {

  public:
    //! Constructor
    MixModel(const int& nscalar, const string& name);

    //! Destructor
    virtual ~MixModel() {}

    //! Interface for echo information on mix model
    virtual void echo() = 0;

  protected:
    const int m_nscalar;           //!< Number of mixing scalars
    const string m_name;           //!< Name of mix model

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
