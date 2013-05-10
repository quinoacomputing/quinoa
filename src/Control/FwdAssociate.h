//******************************************************************************
/*!
  \file      src/Control/FwdAssociate.h
  \author    J. Bakosi
  \date      Thu May  9 19:15:23 2013
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

  // string -> PhysicsType
  using physics_enum = unordered_map< std::string, control::PhysicsType >;
  struct PhysicsEnumStruct {
    static physics_enum make() {
      physics_enum m;
      m["hommix"] = control::PhysicsType::HOMOGENEOUS_MIX;
      m["homhydro"] = control::PhysicsType::HOMOGENEOUS_HYDRO;
      m["spinsflow"] = control::PhysicsType::SPINSFLOW;
      return m;
    }
  };
  static physics_enum PhysicsEnum = PhysicsEnumStruct::make();
//   // ICC: The above struct and definition can be replaced by
//   static physics_enum PhysicsEnum = {
//     { "hommix", control::PhysicsType::HOMOGENEOUS_MIX },
//     { "homhydro", control::PhysicsType::HOMOGENEOUS_HYDRO },
//     { "spinsflow", control::PhysicsType::SPINSFLOW }
//   };

  // string -> PositionType
  using position_enum = unordered_map< std::string, control::PositionType >;
  struct PositionEnumStruct {
    static position_enum make() {
      position_enum m;
      m["invpos"] = control::PositionType::INVISCID;
      m["vispos"] = control::PositionType::VISCOUS;
      return m;
    }
  };
  static position_enum PositionEnum = PositionEnumStruct::make();

  // string -> HydroType
  using hydro_enum = unordered_map< std::string, control::HydroType >;
  struct HydroEnumStruct {
    static hydro_enum make() {
      hydro_enum m;
      m["slm"] = control::HydroType::SLM;
      m["glm"] = control::HydroType::GLM;
      return m;
    }
  };
  static hydro_enum HydroEnum = HydroEnumStruct::make();

  // string -> MixType
  using mix_enum = unordered_map< std::string, control::MixType >;
  struct MixEnumStruct {
    static mix_enum make() {
      mix_enum m;
      m["iem"] = control::MixType::IEM;
      m["iecm"] = control::MixType::IECM;
      m["dir"] = control::MixType::DIRICHLET;
      m["gendir"] = control::MixType::GENERALIZED_DIRICHLET;
      return m;
    }
  };
  static mix_enum MixEnum = MixEnumStruct::make();

} // namespace associate

} // namespace Quinoa

#endif // FwdAssociate_h
