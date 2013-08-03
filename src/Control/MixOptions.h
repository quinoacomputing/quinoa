//******************************************************************************
/*!
  \file      src/Control/MixOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:40:54 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model options and associations
  \details   Mix model options and associations
*/
//******************************************************************************
#ifndef MixOptions_h
#define MixOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Mix model types
enum class MixType : uint8_t { NO_MIX=0,
                               IEM,
                               IECM,
                               DIRICHLET,
                               GENERALIZED_DIRICHLET };

//! Class with base templated on the above enum class with associations
class Mix : public Toggle<MixType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Mix() : Toggle<MixType>(names, values){
      //! Enums -> names
      names[MixType::NO_MIX] = "No mix";
      names[MixType::IEM] = "Interaction by exchange with the mean";
      names[MixType::IECM] =
        "Interaction by exchange with the conditional mean";
      names[MixType::DIRICHLET] = "Dirichlet";
      names[MixType::GENERALIZED_DIRICHLET] = "Generalized Dirichlet";
      //! keywords -> Enums
      values["no_mix"] = MixType::NO_MIX;
      values["mix_iem"] = MixType::IEM;
      values["mix_iecm"] = MixType::IECM;
      values["mix_dir"] = MixType::DIRICHLET;
      values["mix_gendir"] = MixType::GENERALIZED_DIRICHLET;
    }

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;

    std::map<MixType, std::string> names;
    std::map<std::string, MixType> values;
};

} // namespace select

} // namespace Quinoa

#endif // MixOptions_h
