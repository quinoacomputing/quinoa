// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Scheme.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Discretization scheme options for inciter
  \details   Discretization scheme options for inciter
*/
// *****************************************************************************
#ifndef SchemeOptions_h
#define SchemeOptions_h

#include <boost/mpl/vector.hpp>
#include "NoWarning/for_each.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { MatCG
                                , DiagCG
                                , DG };

//! Pack/Unpack SchemeType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, SchemeType& e ) { PUP::pup( p, e ); }

//! \brief Scheme options: outsource to base templated on enum type
class Scheme : public tk::Toggle< SchemeType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::matcg
                                       , kw::diagcg
                                       , kw::dg
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Scheme() :
      tk::Toggle< SchemeType >(
        //! Group, i.e., options, name
        "Discretization scheme",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { SchemeType::MatCG, kw::matcg::name() },
          { SchemeType::DiagCG, kw::diagcg::name() },
          { SchemeType::DG, kw::dg::name() } },
        //! keywords -> Enums
        { { kw::matcg::string(), SchemeType::MatCG },
          { kw::diagcg::string(), SchemeType::DiagCG },
          { kw::dg::string(), SchemeType::DG } } )
    {}

};

} // ctr::
} // inciter::

#endif // SchemeOptions_h
