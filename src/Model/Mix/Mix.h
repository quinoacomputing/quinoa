//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:58:59 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix mode lbase
*/
//******************************************************************************
#ifndef Mix_h
#define Mix_h

#include <string>

#include <QuinoaTypes.h>
#include <Model.h>

namespace Quinoa {

using namespace std;

//! Mix model base
class Mix : public Model {

  public:
    //! Constructor
    Mix(const int& nscalar, const string& name);

    //! Destructor
    virtual ~Mix() {}

    //! Interface for echo information on mix model
    virtual void echo() = 0;

  protected:
    const int m_nscalar;           //!< Number of mixing scalars

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;
};

} // namespace Quinoa

#endif // Mix_h
