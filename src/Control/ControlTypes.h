//******************************************************************************
/*!
  \file      src/Control/ControlTypes.h
  \author    J. Bakosi
  \date      Sat 02 Mar 2013 08:55:18 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for control and parsing
  \details   Types for control and parsing
*/
//******************************************************************************
#ifndef Type_h
#define Type_h

#include <string>
#include <vector>
#include <tuple>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

namespace control {

//! Physics (methods: collection of models) types
enum PhysicsType { NO_PHYSICS=0,
                   HOMOGENEOUS_MIX,
                   SPINSFLOW,
                   NUM_PHYSICS
};

//! Hydrodynamics model types
enum HydroType { NO_HYDRO=0,
                 SLM,
                 GLM,
                 NUM_HYDRO
};

//! Material mix model types
enum MixType { NO_MIX=0,
               IEM,
               IECM,
               DIRICHLET,
               GENERALIZED_DIRICHLET,
               NUM_MIX
};

//! Quantities whose statistics can be estimated
enum Quantity { TRANSPORTED_SCALAR=0,
                VELOCITY,
                PRESSURE,
                DENSITY
};

//! Moment specifies which type of moment is computed for a Quantity in a Term
enum Moment { ORDINARY=0,      //!< Full variable, e.g. Y
              CENTRAL          //!< Fluctuation, e.g. y = Y - <Y>
};

//! Term is a Moment of a Quantity with a field ID to be ensemble averaged
//! The Numbering of field IDs start from 0
//! E.g. 1st pressure fluctuation: {0, PRESSURE, CENTRAL},
//! E.g. mean of 2nd scalar: {1, TRANSPORTED_SCALAR, ORDINARY}
//! readable stores the same information in a more human-readable form
struct Term {
  int field;
  Quantity quantity;
  Moment moment;
  string readable;
};

//! Products are N Terms to be multiplied and ensemble averaged
//! E.g. the scalar flux in x direction needs two terms for ensemble averaging:
//! (Y-<Y>) and (U-<U>), then the moment is <yu> = <(Y-<Y>)(U-<U>)>
//! E.g the third mixed central moment of three scalars needs three terms for
//! ensemble averaging: (Y1-<Y1>), (Y2-<Y2>), and (Y3-<Y3>), then the moment is
//! <y1y2y3> = <(Y1-<Y1>)(Y2-<Y2>)(Y3-<Y3>)>
using Product = vector<Term>;

//! Position enum for accessing fields of tuple Bundle
enum BundlePosition { TITLE=0,
                      PHYSICS,
                      HYDRO,
                      MIX,
                      NSTEP,
                      TERM,
                      DT,
                      NSCALAR,
                      NPAR,
                      TTYI,
                      DUMP,
                      PLTI,
                      PDFI,
                      B,
                      S,
                      KAPPA,
                      C,
                      STATISTICS
};

//! Storage bundle for parsed data
using Bundle = tuple< string,         //!<  0: Title
                      PhysicsType,    //!<  1: Physics
                      HydroType,      //!<  2: Hydrodynamics model
                      MixType,        //!<  3: Material mix model
                      int,            //!<  4: Number of time steps to take
                      real,           //!<  5: Time to terminate time stepping
                      real,           //!<  6: Size of time step
                      int,            //!<  7: Number of mixing scalars
                      int,            //!<  8: Total number of particles
                      int,            //!<  9: TTY output interval
                      int,            //!< 10: Dump output interval
                      int,            //!< 11: Plot output interval
                      int,            //!< 12: PDF output interval
                      vector<real>,   //!< 13: Vector of parameters 'b'
                      vector<real>,   //!< 14: Vector of parameters 'S'
                      vector<real>,   //!< 15: Vector of parameters 'kappa'
                      vector<real>,   //!< 16: Vector of parameters 'c_ij'
                      vector<Product> //!< 17: Vector of statistics
>;

//! Vector of bools indicating whether data is set in Bundle during parsing
using BoolBundle = vector<bool>;

} // namespace control

} // namespace Quinoa

#endif // Type_h
