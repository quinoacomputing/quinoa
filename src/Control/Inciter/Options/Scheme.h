// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Scheme.h
  \copyright 2012-2015, J. Bakosi, 2016-2019, Triad National Security, LLC.
  \brief     Discretization scheme options for inciter
  \details   Discretization scheme options for inciter
*/
// *****************************************************************************
#ifndef SchemeOptions_h
#define SchemeOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"
#include "Centering.h"

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { DiagCG
                                , ALECG
                                , DG
                                , DGP1 
                                , DGP2 };

//! Pack/Unpack SchemeType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, SchemeType& e ) { PUP::pup( p, e ); }

//! \brief Scheme options: outsource to base templated on enum type
class Scheme : public tk::Toggle< SchemeType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::diagcg
                                  , kw::alecg
                                  , kw::dg
                                  , kw::dgp1
                                  , kw::dgp2
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
          { SchemeType::DGP1, kw::dgp1::name() },
          { SchemeType::DGP2, kw::dgp2::name() } },
        //! keywords -> Enums
        { { kw::diagcg::string(), SchemeType::DiagCG },
          { kw::alecg::string(), SchemeType::ALECG },
          { kw::dg::string(), SchemeType::DG },
          { kw::dgp1::string(), SchemeType::DGP1 }, 
          { kw::dgp2::string(), SchemeType::DGP2 } } ) {}

    //! Return scheme centering for SchemeType
    //! \param[in] type Scheme type
    //! \return Mesh centering for scheme type
    tk::Centering centering( SchemeType type ) {
      if ( type == SchemeType::DiagCG ||
           type == SchemeType::ALECG )

        return tk::Centering::NODE;

      else if ( type == SchemeType::DG ||
                type == SchemeType::DGP1 ||
                type == SchemeType::DGP2 )

        return tk::Centering::ELEM;

      else Throw( "No such scheme centering" );
    }
};

} // ctr::
} // inciter::

#endif // SchemeOptions_h
