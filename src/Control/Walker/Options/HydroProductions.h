// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/HydroProductions.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Hydrodynamics production divided by dissipation rate options
  \details   Hydrodynamics production divided by dissipation rate options
*/
// *****************************************************************************
#ifndef HydroProductionsOptions_h
#define HydroProductionsOptions_h

#include <boost/mpl/vector.hpp>
#include "NoWarning/for_each.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"
#include "DiffEq/HydroProductions.h"

namespace walker {
namespace ctr {

//! Hydrodynamics production divided by dissipation rate types
//! \author J. Bakosi
enum class HydroProductionsType : uint8_t { PROD_A005H=0
                                          , PROD_A005S
                                          , PROD_A005L
                                          , PROD_A05H
                                          , PROD_A05S
                                          , PROD_A05L
                                          , PROD_A075H
                                          , PROD_A075S
                                          , PROD_A075L
                                          };

//! \brief Pack/Unpack HydroProductionsType: forward overload to generic enum
//!   class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, HydroProductionsType& e )
{ PUP::pup( p, e ); }

//! HydroProductions options: outsource searches to base templated on enum type
//! \author J. Bakosi
class HydroProductions : public tk::Toggle< HydroProductionsType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::prod_A005H
                                       , kw::prod_A005S
                                       , kw::prod_A005L
                                       , kw::prod_A05H
                                       , kw::prod_A05S
                                       , kw::prod_A05L
                                       , kw::prod_A075H
                                       , kw::prod_A075S
                                       , kw::prod_A075L
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit HydroProductions() :
      tk::Toggle< HydroProductionsType >(
        //! Group, i.e., options, name
        "Hydrodynamics production divided by dissipation rate",
        //! Enums -> names
        { { HydroProductionsType::PROD_A005H, kw::prod_A005H::name() },
          { HydroProductionsType::PROD_A005S, kw::prod_A005S::name() },
          { HydroProductionsType::PROD_A005L, kw::prod_A005L::name() },
          { HydroProductionsType::PROD_A05H, kw::prod_A05H::name() },
          { HydroProductionsType::PROD_A05S, kw::prod_A05S::name() },
          { HydroProductionsType::PROD_A05L, kw::prod_A05L::name() },
          { HydroProductionsType::PROD_A075H, kw::prod_A075H::name() },
          { HydroProductionsType::PROD_A075S, kw::prod_A075S::name() },
          { HydroProductionsType::PROD_A075L, kw::prod_A075L::name() } },
        //! keywords -> Enums
        {  { kw::prod_A005H::string(), HydroProductionsType::PROD_A005H },
           { kw::prod_A005S::string(), HydroProductionsType::PROD_A005S },
           { kw::prod_A005L::string(), HydroProductionsType::PROD_A005L },
           { kw::prod_A05H::string(), HydroProductionsType::PROD_A05H },
           { kw::prod_A05S::string(), HydroProductionsType::PROD_A05S },
           { kw::prod_A05L::string(), HydroProductionsType::PROD_A05L },
           { kw::prod_A075H::string(), HydroProductionsType::PROD_A075H },
           { kw::prod_A075S::string(), HydroProductionsType::PROD_A075S },
           { kw::prod_A075L::string(), HydroProductionsType::PROD_A075L } } ) {}

    //! \brief Return table based on Enum
    //! \param[in] t Enum value of the option requested
    //! \return tk::Table associated to the option
    //! \author J. Bakosi
    tk::Table table( HydroProductionsType t ) const {
      if (t == HydroProductionsType::PROD_A005H)
        return prod_A005H;
      else if (t == HydroProductionsType::PROD_A005S)
        return prod_A005S;
      else if (t == HydroProductionsType::PROD_A005L)
        return prod_A005L;
      else if (t == HydroProductionsType::PROD_A05H)
        return prod_A05H;
      else if (t == HydroProductionsType::PROD_A05S)
        return prod_A05S;
      else if (t == HydroProductionsType::PROD_A05L)
        return prod_A05L;
      else if (t == HydroProductionsType::PROD_A075H)
        return prod_A075H;
      else if (t == HydroProductionsType::PROD_A075S)
        return prod_A075S;
      else if (t == HydroProductionsType::PROD_A075L)
        return prod_A075L;
      else Throw( "Hydrodynamics P/e associated to " +
                  std::to_string( static_cast<uint8_t>(t) ) + " not found" );
    }

};

} // ctr::
} // walker::

#endif // HydroProductionsOptions_h
