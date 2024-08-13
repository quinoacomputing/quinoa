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
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRInitial() :
      tk::Toggle< AMRInitialType >(
        //! Group, i.e., options, name
        "amr_initial",
        //! Enums -> names
        { { AMRInitialType::UNIFORM, "uniform" },
          { AMRInitialType::UNIFORM_DEREFINE, "uniform_derefine" },
          { AMRInitialType::INITIAL_CONDITIONS, "initial_conditions" },
          { AMRInitialType::EDGELIST, "edgelist" },
          { AMRInitialType::COORDINATES, "coords" } },
        //! keywords -> Enums
        { { "uniform", AMRInitialType::UNIFORM },
          { "uniform_derefine", AMRInitialType::UNIFORM_DEREFINE },
          { "initial_conditions", AMRInitialType::INITIAL_CONDITIONS },
          { "edgelist", AMRInitialType::EDGELIST },
          { "coords", AMRInitialType::COORDINATES } } )
    {}
};

} // ctr::
} // inciter::

#endif // InciterAMRInitialOptions_h
