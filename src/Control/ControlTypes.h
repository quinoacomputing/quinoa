//******************************************************************************
/*!
  \file      src/Control/ControlTypes.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 05:22:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for control and parsing
  \details   Types for control and parsing
*/
//******************************************************************************
#ifndef ControlTypes_h
#define ControlTypes_h

#include <string>
#include <vector>
#include <tuple>
#include <ostream>
#include <sstream>

#include <QuinoaTypes.h>

using namespace std;

namespace Quinoa {

namespace control {

//! Physics (methods: collection of models) types
enum PhysicsType { NO_PHYSICS=0,
                   HOMOGENEOUS_MIX,
                   HOMOGENEOUS_HYDRO,
                   SPINSFLOW,
                   NUM_PHYSICS
};

//! Position model types
enum PositionType { NO_POSITION=0,
                    INVISCID,
                    VISCOUS,
                    NUM_POSITION
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
                VELOCITY_X,
                VELOCITY_Y,
                VELOCITY_Z,
                PRESSURE,
                DENSITY
};

//! Moment specifies which type of moment is computed for a Quantity in a Term
enum Moment { ORDINARY=0,      //!< Full variable
              CENTRAL          //!< Fluctuation
};

//! Term is a Moment of a Quantity with a field ID to be ensemble averaged.
//! Internally the Numbering of field IDs starts from 0, but presented to the
//! user as starting from 1. Examples: 1st pressure fluctuation: {0, PRESSURE,
//! CENTRAL, ...}, mean of 2nd scalar: {1, TRANSPORTED_SCALAR, ORDINARY, ...}.
struct Term {
  int field;
  Quantity quantity;
  Moment moment;
  int name;    //! Integer, stores the character code, converted only for output
  bool plot;   //! Shows whether the variable will be plotted
               //! Conceptually, plot should be in Product, since plot will only
               //! be false for a mean that was triggered by a central moment by
               //! one of the Terms of a Product requesting the mean or a model.
               //! However, that would require Product to be a vector<struct>,
               //! which then would need custom comparitors for std::sort() and
               //! std::unique() in Parser::unique(). Since this is not a
               //! performance issue, plot is here in Term.

  //! Constructor
  explicit Term(const int f,
                const Quantity q,
                const Moment m,
                const char n,
                const bool p) : field(f),
                                quantity(q),
                                moment(m),
                                name(n),
                                plot(p) {}

  //! Equal operator for finding unique elements, used by e.g., std::unique()
  //! Test only on field, quantity, and moment
  bool operator== (const Term& term) const {
    if (field == term.field &&
        quantity == term.quantity &&
        moment == term.moment) {
      return true;
    } else {
      return false;
    }
  }

  //! Less than operator for ordering, used by e.g., std::sort().
  //! Test on field, quantity, term, moment, and !plot.
  //! Using operator >, instead of operator <, on plot ensures that if a Term is
  //! user-requested, i.e., plotted, and also triggered by e.g., a model, the
  //! user-requested Term will take precendence.
  bool operator< (const Term& term) const {
    // test on everything except name
    if (field < term.field) {
      return true;
    } else if (field == term.field && quantity < term.quantity) {
      return true;
    } else if (field == term.field && quantity == term.quantity &&
               moment < term.moment) {
      return true;
    } else if (field == term.field && quantity == term.quantity &&
               moment == term.moment && plot > term.plot) {
      return true;
    } else {
      return false;
    }
  }

  //! Operator + for adding Term (name+field ID) to a std::string
  friend string operator+ (const string& lhs, const control::Term& term) {
    stringstream ss;
    ss << lhs << char(term.name) << term.field+1;
    string rhs = ss.str();
    return rhs;
  }

  //! Operator << for writing Term to output streams
  friend ostream& operator<< (ostream& os, const Term& term) {
    os << char(term.name) << term.field+1;
    return os;
  }

  //! Operator << for writing vector<Term> to output streams
  friend ostream& operator<< (ostream& os, const vector<Term>& vec) {
    os << " <";
    for (auto& w : vec) os << w;
    os << ">";
    return os;
  }

  //! Operator <<= for writing requested vector<Term> to output streams
  friend ostream& operator<<= (ostream& os, const vector<Term>& vec) {
    if (vec[0].plot) {
      os << " <";
      for (auto& w : vec) os << w;
      os << ">";
    }
    return os;
  }
};

//! Products are N Terms to be multiplied and ensemble averaged
//! E.g. the scalar flux in x direction needs two terms for ensemble averaging:
//! (Y-<Y>) and (U-<U>), then the moment is <yu> = <(Y-<Y>)(U-<U>)>
//! E.g the third mixed central moment of three scalars needs three terms for
//! ensemble averaging: (Y1-<Y1>), (Y2-<Y2>), and (Y3-<Y3>), then the moment is
//! <y1y2y3> = <(Y1-<Y1>)(Y2-<Y2>)(Y3-<Y3>)>
using Product = vector<Term>;

//! Position enum for accessing fields of tuple Bundle using names as in struct
enum BundlePosition { TITLE=0,
                      PHYSICS,
                      POSITION,
                      HYDRO,
                      MIX,
                      NSTEP,
                      TERM,
                      DT,
                      NPOSITION,
                      NVELOCITY,
                      NSCALAR,
                      NPAR,
                      TTYI,
                      DUMP,
                      PLTI,
                      PDFI,
                      GLOB,
                      PDFNAME,
                      GLOBNAME,
                      PLOTNAME,
                      B,
                      S,
                      KAPPA,
                      C,
                      C0,
                      STATISTICS
};

//! Storage bundle for parsed data
using Bundle = tuple<
  string,               //!< Problem Title
  PhysicsType,          //!< Selected physics
  PositionType,         //!< Selected position model
  HydroType,            //!< Selected hydrodynamics model
  MixType,              //!< Selected material mix model
  int,                  //!< Number of time steps to take
  real,                 //!< Time to terminate time stepping
  real,                 //!< Size of time step
  int,                  //!< Number of position components in position model
  int,                  //!< Number of velocity components in hydro model
  int,                  //!< Number of mixing scalars in material mix model
  int,                  //!< Total number of particles
  int,                  //!< TTY output interval
  int,                  //!< Dump output interval
  int,                  //!< Plot output interval
  int,                  //!< PDF output interval
  int,                  //!< Glob output interval
  string,               //!< PDF base filename
  string,               //!< Glob filename
  string,               //!< Plot base filename
  vector<real>,         //!< Parameters 'b' in Dirichlet mix models
  vector<real>,         //!< Parameters 'S' in Dirichlet mix models
  vector<real>,         //!< Parameters 'kappa' in Dirichlet mix models
  vector<real>,         //!< Parameters 'c_ij' in GenDirichlet mix models
  real,                 //!< Parameter C0 in the simplified Langevin model
  vector<Product>       //!< Requested (and triggered) statistics
>;

//! Vector of bools indicating whether data is set in Bundle during parsing
using BoolBundle = vector<bool>;

} // namespace control

} // namespace Quinoa

#endif // ControlTypes_h
