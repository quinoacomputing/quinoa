//******************************************************************************
/*!
  \file      src/Model/MixModel/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Thu Nov 15 15:50:11 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <MixModel.h>

namespace Quinoa {

//! GeneralizedDirichlet : MixModel
class GeneralizedDirichlet : public MixModel {

  public:
    //! Constructor
    GeneralizedDirichlet(const int& nscalar);

    //! Destructor
    virtual ~GeneralizedDirichlet() {}

    //! Echo information on Generalized Dirichlet model
    virtual void echo();

    //! Initialize Generalized Dirichlet model
    virtual void init();

    //! Set initial conditions
    void setIC();

  private:
    //! Don't permit copy constructor
    GeneralizedDirichlet(const GeneralizedDirichlet&) = delete;
    //! Don't permit copy assigment
    GeneralizedDirichlet& operator=(const GeneralizedDirichlet&) = delete;
    //! Don't permit move constructor
    GeneralizedDirichlet(GeneralizedDirichlet&&) = delete;
    //! Don't permit move assigment
    GeneralizedDirichlet& operator=(GeneralizedDirichlet&&) = delete;
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
