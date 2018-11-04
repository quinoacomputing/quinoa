// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Scheme.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { MatCG
                                , DiagCG
                                , DG
                                , DGP1 };

//! Scheme centering types
enum class Centering : uint8_t { NODE
                               , ELEM };

//! Pack/Unpack SchemeType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, SchemeType& e ) { PUP::pup( p, e ); }

//! \brief Scheme options: outsource to base templated on enum type
class Scheme : public tk::Toggle< SchemeType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::matcg
                                  , kw::diagcg
                                  , kw::dg
                                  , kw::dgp1
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Scheme() :
      tk::Toggle< SchemeType >(
        //! Group, i.e., options, name
        kw::scheme::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { SchemeType::MatCG, kw::matcg::name() },
          { SchemeType::DiagCG, kw::diagcg::name() },
          { SchemeType::DG, kw::dg::name() },
          { SchemeType::DGP1, kw::dgp1::name() } },
        //! keywords -> Enums
        { { kw::matcg::string(), SchemeType::MatCG },
          { kw::diagcg::string(), SchemeType::DiagCG },
          { kw::dg::string(), SchemeType::DG },
          { kw::dgp1::string(), SchemeType::DGP1 } } ) {}

    //! Return scheme centering for SchemeType
    //! \param[in] type Scheme type
    //! \return Centering for scheme type
    Centering centering( SchemeType type ) {
      if (type == SchemeType::MatCG || type == SchemeType::DiagCG)
        return Centering::NODE;
      else if (type == SchemeType::DG || type == SchemeType::DGP1)
        return Centering::ELEM;
      else Throw( "No such scheme centering" );
    }
};

} // ctr::
} // inciter::

#endif // SchemeOptions_h
