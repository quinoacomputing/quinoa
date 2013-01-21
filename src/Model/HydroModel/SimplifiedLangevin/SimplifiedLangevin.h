//******************************************************************************
/*!
  \file      src/Model/HydroModel/SimplifiedLangevin/SimplifiedLangevin.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:48:21 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SimplifiedLangevin_h
#define SimplifiedLangevin_h

#include <QuinoaTypes.h>
#include <HydroModel.h>

namespace Quinoa {

//! SimplifiedLangevin : HydroModel
class SimplifiedLangevin : public HydroModel {

  public:
    //! Constructor
    SimplifiedLangevin();

    //! Destructor
    virtual ~SimplifiedLangevin() {}

    //! Echo information on the simplified Langevin model
    virtual void echo();

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;
};

} // namespace Quinoa

#endif // SimplifiedLangevin_h
