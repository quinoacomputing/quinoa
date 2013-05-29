//******************************************************************************
/*!
  \file      src/Control/MixOptions.h
  \author    J. Bakosi
  \date      Wed May 29 07:34:24 2013
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
enum class MixTypes : uint8_t { NO_MIX=0,
                                IEM,
                                IECM,
                                DIRICHLET,
                                GENERALIZED_DIRICHLET };

//! Class with base templated on the above enum class with associations
class Mix : public Toggle<MixTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Mix() : Toggle<MixTypes>(names, values){
      //! Enums -> names
      names[MixTypes::NO_MIX] = "No mix";
      names[MixTypes::IEM] = "Interaction by exchange with the mean";
      names[MixTypes::IECM] =
        "Interaction by exchange with the conditional mean";
      names[MixTypes::DIRICHLET] = "Dirichlet";
      names[MixTypes::GENERALIZED_DIRICHLET] = "Dirichlet";
      //! keywords -> Enums
      values["no_mix"] = MixTypes::NO_MIX;
      values["iem"] = MixTypes::IEM;
      values["iecm"] = MixTypes::IECM;
      values["dir"] = MixTypes::DIRICHLET;
      values["gendir"] = MixTypes::GENERALIZED_DIRICHLET;
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

    std::map<MixTypes, std::string> names;
    std::map<std::string, MixTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // MixOptions_h
