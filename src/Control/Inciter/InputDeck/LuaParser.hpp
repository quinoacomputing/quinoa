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

#include "Inciter/CmdLine/CmdLine.hpp"
#include "InputDeck.hpp"

namespace tk { class Print; }

namespace inciter {

//! \brief Control file lua-parser for Inciter.
//! \details This class is used to interface with sol2, for the purpose of
//!   parsing the control file for the computational shock hydrodynamics tool,
//!   Inciter.
class LuaParser {

  public:
    //! Constructor
    explicit LuaParser( const tk::Print& print,
                              const ctr::CmdLine& cmdline,
                              ctr::InputDeck& inputdeck );

    //! Store lua inputdeck in custom struct
    void storeInputDeck(
      const sol::table& lua_ideck,
      ctr::InputDeck& gideck );

    //! Check and store material property into inpudeck storage
    void checkStoreMatProp(
      const sol::table table,
      const std::string key,
      std::size_t vecsize,
      std::vector< tk::real >& storage );

    //! Check and store material property vector into inpudeck storage
    void checkStoreMatPropVec(
      const sol::table table,
      const std::string key,
      std::size_t nspec,
      std::size_t vecsize,
      std::vector<std::vector< tk::real >>& storage );

    //! Check and store material property vector of vectors into inpudeck storage
    void checkStoreMatPropVecVec(
      const sol::table table,
      const std::string key,
      std::size_t nspec,
      std::size_t vecsize1,
      std::size_t vecsize2,
      std::vector<std::vector<std::vector< tk::real >>>& storage );

    //! Check and store field output variables
    void addOutVar(
      const std::string& varname,
      const std::string& alias,
      std::vector< char >& depv,
      std::size_t nmat,
      std::size_t nspec,
      inciter::ctr::PDEType pde,
      tk::Centering c,
      std::vector< inciter::ctr::OutVar >& foutvar );

  private:
    const std::string m_filename;             //!< Name of file to parse

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

    //! \brief Check validity of keywords within a particular block
    //! \tparam tags Tags addressing the said block in the input deck
    //! \param[in] block Sol table of the input deck block read in from the
    //!   lua file
    //! \param[in] blk_name Name of the block, for error clarity
    template< typename... tags >
    void checkBlock(const sol::table& block, const std::string& blk_name) const
    {
      for (const auto& kvp : block) {
        bool is_valid(false);
        auto ukw = kvp.first.as<std::string>();
        brigand::for_each< tags... >( checkKw(ukw, is_valid) );
        if (!is_valid)
          Throw("Invalid keyword '" + ukw + "' in '" + blk_name + "' block.");
      }
    }

    // Check if a keyword matches the existing ones
    struct checkKw {
      std::string user_kw;
      // reference to bool keeping track of kw-match
      bool& is_valid;
      // Constructor
      checkKw( const std::string& ukw, bool& isv ) :
        user_kw(ukw), is_valid(isv) {}
      //! Function to call for each keyword type
      template< typename U > void operator()( brigand::type_<U> ) {
        auto spec_key = U::name();
        // only check if not previously matched
        if (!is_valid) {
          if (user_kw == spec_key) is_valid = true;
          else is_valid = false;
        }
      }
    };
};

} // namespace inciter

#endif // InciterLuaParser_h
