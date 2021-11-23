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
    {}
};

} // ctr::
} // inciter::

#endif // InciterAMRInitialOptions_h
