//******************************************************************************
/*!
  \file      src/Control/FwdAssociate.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 12:49:18 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Forward, keyword (string) to data associations
  \details   Forward, keyword (string) to data associations
*/
//******************************************************************************
#ifndef FwdAssociate_h
#define FwdAssociate_h

#include <unordered_map>

namespace Quinoa {

namespace associate {

  // Editing anything below should be accompanied by the corresponding changes
  // in BackAssociate.h as well.

  using namespace control;

  // string -> PhysicsType
  using physics_enum = unordered_map< std::string, PhysicsType >;
  struct PhysicsEnumStruct {
    static physics_enum make() {
      physics_enum m;
      m["homdir"] = PhysicsType::HOMOGENEOUS_DIRICHLET;
      m["homgendir"] = PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET;
      m["spinsflow"] = PhysicsType::SPINSFLOW;
      return m;
    }
  };
  static physics_enum PhysicsEnum = PhysicsEnumStruct::make();
//   // ICC: The above struct and definition can be replaced by
//   static physics_enum PhysicsEnum = {
//     { "homdir", PhysicsType::HOMOGENEOUS_DIRICHLET },
//     { "homgendir", PhysicsType::HOMOGENEOUS_GENERALIZED_DIRICHLET },
//     { "spinsflow", PhysicsType::SPINSFLOW }
//   };

  // string -> HydroType
  using hydro_enum = unordered_map< std::string, HydroType >;
  struct HydroEnumStruct {
    static hydro_enum make() {
      hydro_enum m;
      m["slm"] = HydroType::SLM;
      m["glm"] = HydroType::GLM;
      return m;
    }
  };
  static hydro_enum HydroEnum = HydroEnumStruct::make();

  // string -> MixType
  using mix_enum = unordered_map< std::string, MixType >;
  struct MixEnumStruct {
    static mix_enum make() {
      mix_enum m;
      m["iem"] = MixType::IEM;
      m["iecm"] = MixType::IECM;
      m["dir"] = MixType::DIRICHLET;
      m["gendir"] = MixType::GENERALIZED_DIRICHLET;
      return m;
    }
  };
  static mix_enum MixEnum = MixEnumStruct::make();

} // namespace associate

} // namespace Quinoa

#endif // FwdAssociate_h
