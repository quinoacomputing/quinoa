//******************************************************************************
/*!
  \file      src/Parser/Keywords.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 08:09:32 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Keywords
  \details   All keywords recognized by the parser
*/
//******************************************************************************
#ifndef Keywords_h
#define Keywords_h

// Keywords accepted by the parser
namespace keyword {

  // Problem title
  using title = pegtl::string<t,i,t,l,e>;

  // End of block
  using end = pegtl::string< e,n,d >;

  // Select physics:
  //   * Homogeneous Dirichlet
  using homdir = pegtl::string< h,o,m,d,i,r >;
  //   * Homogeneous generalized Dirichlet
  using homgendir = pegtl::string< h,o,m,g,e,n,d,i,r >;
  //   * Standalone-particle incompressible Navier-Stokes flow
  using spinsflow = pegtl::string< s,p,i,n,s,f,l,o,w >;

  // Select hydrodynamics model:
  //   * Simplified Langevin model
  using slm = pegtl::string< s,l,m >;
  //   * Generalized Langevin model
  using glm = pegtl::string< g,l,m >;

  // Select material mix model:
  //   * Interaction by exchange with the mean
  using iem = pegtl::string< i,e,m >;
  //   * Interaction by exchange with the conditional mean
  using iecm = pegtl::string< i,e,c,m >;
  //   * Dirichlet
  using dir = pegtl::string< d,i,r >;
  //   * generalized Dirichlet
  using gendir = pegtl::string< g,e,n,d,i,r >;

  // Number of time steps to take
  using nstep = pegtl::string< n,s,t,e,p >;

  // Terminate time stepping at this value
  using term = pegtl::string< t, e, r, m >;

  // Size of time step
  using dt = pegtl::string< d,t >;

  // Start of hydrodynamics model specification block
  using hydro = pegtl::string< h,y,d,r,o >;

  // Start of material mix model specification block
  using mix = pegtl::string< m,i,x >;

  // Number of mixing scalars
  using nscalar = pegtl::string< n,s,c,a,l,a,r >;

  // Total number of particles
  using npar = pegtl::string< n,p,a,r >;

  // One-line screen output every few time step
  using echo = pegtl::string< e,c,h,o >;

} // namespace keyword


// Associations between parsed strings and Quinoa enums
namespace associate {

  using namespace control;

  struct PhysicsMap {
    static unordered_map< std::string, PhysicsType > make() {
      unordered_map< std::string, PhysicsType > m;
      m["homdir"] = PhysicsType::HOMOGENEOUS_DIRICHLET;
      m["homgendir"] = PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET;
      m["spinsflow"] = PhysicsType::SPINSFLOW;
      return m;
    }
  };
  unordered_map< std::string, PhysicsType > Physics = PhysicsMap::make();
// ICC: The above can be replaced by the below
//   unordered_map<std::string, PhysicsType> Physics = {
//     { "homdir", PhysicsType::HOMOGENEOUS_DIRICHLET },
//     { "homgendir", PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET }
//     { "spinsflow", PhysicsType::SPINSFLOW }
//   };

  struct HydroMap {
    static unordered_map< std::string, HydroType > make() {
      unordered_map< std::string, HydroType > m;
      m["slm"] = HydroType::SLM;
      m["glm"] = HydroType::GLM;
      return m;
    }
  };
  unordered_map< std::string, HydroType > Hydro = HydroMap::make();
// ICC: The above can be replaced by the below
//   unordered_map<std::string, HydroType> Hydro = {
//     { "slm", HydroType::SLM },
//     { "glm", HydroType::GLM }
//   };

  struct MixMap {
    static unordered_map< std::string, MixType > make() {
      unordered_map< std::string, MixType > m;
      m["iem"]    = MixType::IEM;
      m["iecm"]   = MixType::IECM;
      m["dir"]    = MixType::DIRICHLET;
      m["gendir"] = MixType::GENERALIZED_DIRICHLET;
      return m;
    }
  };
  unordered_map< std::string, MixType > Mix = MixMap::make();
// ICC: The above can be replaced by the below
//   unordered_map<std::string, MixType> Mix = {
//     { "iem",       MixType::IEM },
//     { "iecm",      MixType::IECM },
//     { "dir",       MixType::DIRICHLET },
//     { "gendirdir", MixType::GENERALIZED_DIRICHLET }
//   };

} // namespace associate

#endif // Keywords_h
