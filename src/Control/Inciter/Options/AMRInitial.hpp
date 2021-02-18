// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/AMRInitial.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Initial (before t=0) adaptive mesh refinement (AMR) options
  \details   Initial (before t=0) adaptive mesh refinement (AMR) options
*/
// *****************************************************************************
#ifndef InciterAMRInitialOptions_h
#define InciterAMRInitialOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Initial AMR types
enum class AMRInitialType : uint8_t { UNIFORM
                                    , UNIFORM_DEREFINE
                                    , INITIAL_CONDITIONS
                                    , EDGELIST
                                    , COORDINATES };

//! Pack/Unpack AMRInitialType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, AMRInitialType& e )
{ PUP::pup( p, e ); }

//! AMRInitial options: outsource searches to base templated on enum type
class AMRInitial : public tk::Toggle< AMRInitialType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::amr_uniform
                                  , kw::amr_uniform_derefine
                                  , kw::amr_initial_conditions
                                  , kw::amr_edgelist
                                  , kw::amr_coords >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRInitial() :
      tk::Toggle< AMRInitialType >(
        //! Group, i.e., options, name
        kw::amr_initial::name(),
        //! Enums -> names
        { { AMRInitialType::UNIFORM, kw::amr_uniform::name() },
          { AMRInitialType::UNIFORM_DEREFINE,
            kw::amr_uniform_derefine::name() },
          { AMRInitialType::INITIAL_CONDITIONS,
            kw::amr_initial_conditions::name() },
          { AMRInitialType::EDGELIST, kw::amr_edgelist::name() },
          { AMRInitialType::COORDINATES, kw::amr_coords::name() } },
        //! keywords -> Enums
        { { kw::amr_uniform::string(), AMRInitialType::UNIFORM },
          { kw::amr_uniform_derefine::string(),
            AMRInitialType::UNIFORM_DEREFINE },
          { kw::amr_initial_conditions::string(),
            AMRInitialType::INITIAL_CONDITIONS },
          { kw::amr_edgelist::string(), AMRInitialType::EDGELIST },
          { kw::amr_coords::string(), AMRInitialType::COORDINATES } } )
    {
       brigand::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the option requested
    //! \return Policy code of the option
    const std::string& code( AMRInitialType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for initial AMR type \"")
              << p << "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( brigand::type_<U> ) {
        static_assert( tk::HasTypedef_code_v< typename U::info >,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< AMRInitialType, std::string > policy {
        { AMRInitialType::UNIFORM, *kw::amr_uniform::code() }
      , { AMRInitialType::UNIFORM_DEREFINE, *kw::amr_uniform_derefine::code() }
      , { AMRInitialType::INITIAL_CONDITIONS,
          *kw::amr_initial_conditions::code() }
      , { AMRInitialType::EDGELIST, *kw::amr_edgelist::code() }
      , { AMRInitialType::COORDINATES, *kw::amr_coords::code() }
    };
};

} // ctr::
} // inciter::

#endif // InciterAMRInitialOptions_h
