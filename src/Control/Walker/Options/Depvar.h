// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/Depvar.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Dependent variable options for walker
  \details   Dependent variable options for walker
*/
// *****************************************************************************
#ifndef DepvarOptions_h
#define DepvarOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace walker {
namespace ctr {

//! Dependent variable options types
enum class DepvarType : uint8_t { FULLVAR=0
                                , FLUCTUATION };

//! Pack/Unpack DepvarType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, DepvarType& e ) { PUP::pup( p, e ); }

//! Dependent variable options: outsource to base templated on enum type
class Depvar : public tk::Toggle< DepvarType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::fullvar
                                  , kw::fluctuation
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Depvar() :
      tk::Toggle< DepvarType >(
        //! Group, i.e., options, name
        "Solver for",
        //! Enums -> names
        { { DepvarType::FULLVAR, kw::fullvar::name() },
          { DepvarType::FLUCTUATION, kw::fluctuation::name() } },
        //! keywords -> Enums
        { { kw::fullvar::string(), DepvarType::FULLVAR },
          { kw::fluctuation::string(), DepvarType::FLUCTUATION } } )
    {}
};

} // ctr::
} // walker::

#endif // DepvarOptions_h
