// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/Depvar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Dependent variable options for walker
  \details   Dependent variable options for walker
*/
// *****************************************************************************
#ifndef DepvarOptions_h
#define DepvarOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace walker {
namespace ctr {

//! Dependent variable options types
enum class DepvarType : uint8_t { FULLVAR=0
                                , FLUCTUATION
                                , PRODUCT
                                , FLUCTUATING_MOMENTUM };

//! Pack/Unpack DepvarType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, DepvarType& e ) { PUP::pup( p, e ); }

//! Dependent variable options: outsource to base templated on enum type
class Depvar : public tk::Toggle< DepvarType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::fullvar
                                  , kw::fluctuation
                                  , kw::product
                                  , kw::fluctuating_momentum
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Depvar() :
      tk::Toggle< DepvarType >(
        //! Group, i.e., options, name
        "solve for",
        //! Enums -> names
        { { DepvarType::FULLVAR, kw::fullvar::name() },
          { DepvarType::FLUCTUATION, kw::fluctuation::name() },
          { DepvarType::PRODUCT, kw::product::name() },
          { DepvarType::FLUCTUATING_MOMENTUM, kw::fluctuating_momentum::name() }
        },
        //! keywords -> Enums
        { { kw::fullvar::string(), DepvarType::FULLVAR },
          { kw::fluctuation::string(), DepvarType::FLUCTUATION },
          { kw::product::string(), DepvarType::PRODUCT },
          { kw::fluctuating_momentum::string(),
            DepvarType::FLUCTUATING_MOMENTUM }
        } )
    {}
};

} // ctr::
} // walker::

#endif // DepvarOptions_h
