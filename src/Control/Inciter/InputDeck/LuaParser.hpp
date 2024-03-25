// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/LuaParser.hpp
  \copyright 2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's lua input deck file parser
  \details   This file declares the input deck, i.e., control file, parser for
    the computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InciterLuaParser_h
#define InciterLuaParser_h

#include "NoWarning/sol.hpp"

#include "FileParser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "New2InputDeck.hpp"

namespace tk { class Print; }

namespace inciter {

//! \brief Control file lua-parser for Inciter.
//! \details This class is used to interface with sol2, for the purpose of
//!   parsing the control file for the computational shock hydrodynamics tool,
//!   Inciter.
class LuaParser : public tk::FileParser {

  public:
    //! Constructor
    explicit LuaParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::New2InputDeck& inputdeck );

    //! Store lua inputdeck in custom struct
    void storeInputDeck(
      const sol::table& lua_ideck,
      ctr::New2InputDeck& gideck );

    //! Check and store material property into inpudeck storage
    void checkStoreMatProp(
      const sol::table table,
      const std::string key,
      std::size_t vecsize,
      std::vector< tk::real >& storage );

    //! Check and store field output variables
    void addOutVar(
      const std::string& varname,
      std::vector< char >& depv,
      std::size_t nmat,
      inciter::ctr::PDEType pde,
      tk::Centering c,
      std::vector< inciter::ctr::OutVar >& foutvar );

  private:
    //! Assign parameter to inputdeck entry if specified, else default
    //! \tparam N Type of parameter being read/assigned
    //! \param[in] table Sol-table which contains said parameter
    //! \param[in] key Key for said parameter in Sol-table
    //! \param[in,out] storage Storage space in inputdeck where said parameter
    //!   is to be stored
    //! \param[in] dflt Default value of said parameter, if unspecified
    template< typename N > void
    storeIfSpecd(
      const sol::table table,
      const std::string key,
      N& storage,
      const N dflt )
    {
      auto sol_var = table[key];
      if (sol_var.valid())
        storage = sol_var;
      else
        storage = dflt;
    }

    //! Assign Option to inputdeck entry if specified, else default
    //! \tparam Op Option-Type of parameter being read/assigned
    //! \tparam OpClass Option-class of parameter being read/assigned
    //! \param[in] table Sol-table which contains said parameter
    //! \param[in] key Key for said parameter in Sol-table
    //! \param[in,out] storage Storage space in inputdeck where said parameter
    //!   is to be stored
    //! \param[in] dflt Default value of said parameter, if unspecified
    template< typename Op, class OpClass > void
    storeOptIfSpecd(
      const sol::table table,
      const std::string key,
      Op& storage,
      const Op dflt )
    {
      OpClass opt;
      auto sol_var = table[key];
      if (sol_var.valid())
        storage = opt.value(sol_var);
      else
        storage = dflt;
    }

    //! Assign vector parameter to inputdeck entry if specified, else default
    //! \tparam N Type of parameter vector being read/assigned
    //! \param[in] table Sol-table which contains said parameter
    //! \param[in] key Key for said parameter in Sol-table
    //! \param[in,out] storage Storage space in inputdeck where said parameter
    //!   is to be stored
    //! \param[in] dflt Default value of said parameter, if unspecified
    template< typename N > void
    storeVecIfSpecd(
      const sol::table table,
      const std::string key,
      std::vector< N >& storage,
      const std::vector< N >& dflt )
    {
      auto sol_vec = table[key];
      if (sol_vec.valid()) {
        for (std::size_t i=0; i<sol::table(sol_vec).size(); ++i)
          storage.push_back(sol_vec[i+1]);
      }
      else
        storage = dflt;
    }

    //! Assign vector of Options to inputdeck entry if specified, else default
    //! \tparam Op Option-Type of parameter being read/assigned
    //! \tparam OpClass Option-class of parameter being read/assigned
    //! \param[in] table Sol-table which contains said parameter
    //! \param[in] key Key for said parameter in Sol-table
    //! \param[in,out] storage Storage space in inputdeck where said parameter
    //!   is to be stored
    //! \param[in] dflt Default value of said parameter, if unspecified
    template< typename Op, class OpClass > void
    storeOptVecIfSpecd(
      const sol::table table,
      const std::string key,
      std::vector< Op >& storage,
      const std::vector< Op >& dflt )
    {
      OpClass opt;
      auto sol_vec = table[key];
      if (sol_vec.valid()) {
        for (std::size_t i=0; i<sol::table(sol_vec).size(); ++i)
          storage.push_back(opt.value(sol_vec[i+1]));
      }
      else
        storage = dflt;
    }
};

} // namespace inciter

#endif // InciterLuaParser_h
