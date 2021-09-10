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
#include "Keywords.hpp"
#include "PUPUtil.hpp"
#include "Centering.hpp"

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { DiagCG
                                , ALECG
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
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::diagcg
                                  , kw::alecg
                                  , kw::dg
                                  , kw::p0p1
                                  , kw::dgp1
                                  , kw::dgp2
                                  , kw::pdg
                                  , kw::fv
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Scheme() :
      tk::Toggle< SchemeType >(
        //! Group, i.e., options, name
        kw::scheme::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { SchemeType::DiagCG, kw::diagcg::name() },
          { SchemeType::ALECG, kw::alecg::name() },
          { SchemeType::DG, kw::dg::name() },
          { SchemeType::P0P1, kw::p0p1::name() },
          { SchemeType::DGP1, kw::dgp1::name() },
          { SchemeType::DGP2, kw::dgp2::name() },
          { SchemeType::PDG, kw::pdg::name() },
          { SchemeType::FV, kw::fv::name() } },
        //! keywords -> Enums
        { { kw::diagcg::string(), SchemeType::DiagCG },
          { kw::alecg::string(), SchemeType::ALECG },
          { kw::dg::string(), SchemeType::DG },
          { kw::p0p1::string(), SchemeType::P0P1 }, 
          { kw::dgp1::string(), SchemeType::DGP1 }, 
          { kw::dgp2::string(), SchemeType::DGP2 },
          { kw::pdg::string(), SchemeType::PDG },
          { kw::fv::string(), SchemeType::FV } } ) {}

    //! Return scheme centering for SchemeType
    //! \param[in] type Scheme type
    //! \return Mesh centering for scheme type
    tk::Centering centering( SchemeType type ) {
      if ( type == SchemeType::DiagCG ||
           type == SchemeType::ALECG )

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
