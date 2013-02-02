//******************************************************************************
/*!
  \file      src/Control/BackAssociate.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 12:49:31 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Backward, data to name (string) associations
  \details   Backward, data to name (string) associations
*/
//******************************************************************************
#ifndef BackAssociate_h
#define BackAssociate_h

#include <map>

namespace Quinoa {

namespace associate {

  // Editing anything below should be accompanied by the corresponding changes
  // in FwdAssociate.h as well.

  using namespace control;

  // PhysicsType -> string
  using physics_name = map< PhysicsType, std::string >;
  struct PhysicsNameStruct {
    static physics_name make() {
      physics_name m;
      m[PhysicsType::HOMOGENEOUS_DIRICHLET] = "homdir";
      m[PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET] = "homgendir";
      m[PhysicsType::SPINSFLOW] = "spinsflow";
      return m;
    }
  };
  static physics_name PhysicsName = PhysicsNameStruct::make();
//   // ICC: The above struct and definition can be replaced by
//   static physics_name PhysicsName = {
//     { PhysicsType::HOMOGENEOUS_DIRICHLET, "homdir" },
//     { PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET, "homgendir" },
//     { PhysicsType::SPINSFLOW, "spinsflow" }
//   };

  // HydroType -> string
  using hydro_name = map< HydroType, std::string >;
  struct HydroNameStruct {
    static hydro_name make() {
      hydro_name m;
      m[HydroType::SLM] = "slm";
      m[HydroType::GLM] = "glm";
      return m;
    }
  };
  static hydro_name HydroName = HydroNameStruct::make();

  // MixType -> string
  using mix_name = map< MixType, std::string >;
  struct MixNameStruct {
    static mix_name make() {
      mix_name m;
      m[MixType::IEM] = "iem";
      m[MixType::IECM] = "iecm";
      m[MixType::DIRICHLET] = "dir";
      m[MixType::GENERALIZED_DIRICHLET] = "gendir";
      return m;
    }
  };
  static mix_name MixName = MixNameStruct::make();

} // namespace associate

} // namespace Quinoa

#endif // BackAssociate_h
