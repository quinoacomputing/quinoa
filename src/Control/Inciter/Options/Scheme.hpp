// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Scheme.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Discretization scheme options for inciter
  \details   Discretization scheme options for inciter
*/
// *****************************************************************************
#ifndef SchemeOptions_h
#define SchemeOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"
#include "Centering.hpp"

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { ALECG
                                , OversetFE
                                , DG
                                , P0P1 
                                , DGP1 
                                , DGP2
                                , PDG
                                , FV };

//! Pack/Unpack SchemeType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, SchemeType& e ) { PUP::pup( p, e ); }

//! \brief Scheme options: outsource to base templated on enum type
class Scheme : public tk::Toggle< SchemeType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Scheme() :
      tk::Toggle< SchemeType >(
        //! Group, i.e., options, name
        "Scheme",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { SchemeType::ALECG, "alecg" },
          { SchemeType::OversetFE, "oversetfe" },
          { SchemeType::DG, "dg" },
          { SchemeType::P0P1, "p0p1" },
          { SchemeType::DGP1, "dgp1" },
          { SchemeType::DGP2, "dgp2" },
          { SchemeType::PDG, "pdg" },
          { SchemeType::FV, "fv" } },
        //! keywords -> Enums
        { { "alecg", SchemeType::ALECG },
          { "oversetfe", SchemeType::OversetFE },
          { "dg", SchemeType::DG },
          { "p0p1", SchemeType::P0P1 }, 
          { "dgp1", SchemeType::DGP1 }, 
          { "dgp2", SchemeType::DGP2 },
          { "pdg", SchemeType::PDG },
          { "fv", SchemeType::FV } } ) {}

    //! Return scheme centering for SchemeType
    //! \param[in] type Scheme type
    //! \return Mesh centering for scheme type
    tk::Centering centering( SchemeType type ) {
      if ( type == SchemeType::ALECG ||
           type == SchemeType::OversetFE )

        return tk::Centering::NODE;

      else if ( type == SchemeType::DG ||
                type == SchemeType::P0P1 ||
                type == SchemeType::DGP1 ||
                type == SchemeType::DGP2 ||
                type == SchemeType::PDG ||
                type == SchemeType::FV )

        return tk::Centering::ELEM;

      else Throw( "No such scheme centering" );
    }
};

} // ctr::
} // inciter::

#endif // SchemeOptions_h
