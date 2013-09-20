//******************************************************************************
/*!
  \file      src/Control/MixOptions.h
  \author    J. Bakosi
  \date      Fri Sep 20 13:37:56 2013
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
#include <QuinoaKeywords.h>

namespace quinoa {
namespace sel {

//! Mix model types
enum class MixType : uint8_t { NO_MIX=0,
                               IEM,
                               IECM,
                               DIRICHLET,
                               GENERALIZED_DIRICHLET };

//! Class with base templated on the above enum class with associations
class Mix : public Toggle<MixType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Mix() : Toggle<MixType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;

    //! Enums -> names
    const std::map<MixType, std::string> names {
      { MixType::NO_MIX, "" },
      { MixType::IEM, "Interaction by exchange with the mean" },
      { MixType::IECM, "Interaction by exchange with the conditional mean" },
      { MixType::DIRICHLET, "Dirichlet" },
      { MixType::GENERALIZED_DIRICHLET, "Generalized Dirichlet" }
    };

    //! Get access to mix keywords
    const grm::kw::mix_iem iem {};
    const grm::kw::mix_iecm iecm {};
    const grm::kw::mix_dir dir {};
    const grm::kw::mix_gendir gendir {};

    //! keywords -> Enums
    const std::map<std::string, MixType> values {
      { "no_mix", MixType::NO_MIX },
      { iem.string(), MixType::IEM },
      { iecm.string(), MixType::IECM },
      { dir.string(), MixType::DIRICHLET },
      { gendir.string(), MixType::GENERALIZED_DIRICHLET }
    };
};

} // sel::
} // quinoa::

#endif // MixOptions_h
