// *****************************************************************************
/*!
  \file      src/Control/Options/PartitioningAlgorithm.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Mesh partitioning algorithm options
  \details   Mesh partitioning algorithm options
*/
// *****************************************************************************
#ifndef InciterPartitioningAlgorithmOptions_h
#define InciterPartitioningAlgorithmOptions_h

#include <map>

#include "NoWarning/vector.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! Mesh partitioning algorithm types
enum class PartitioningAlgorithmType : uint8_t { RCB,
                                                 RIB,
                                                 HSFC,
                                                 MJ,
                                                 PHG };

//! \brief Pack/Unpack PartitioningAlgorithmType: forward overload to generic
//!   enum class packer
inline void operator|( PUP::er& p, PartitioningAlgorithmType& e )
{ PUP::pup( p, e ); }

//! \brief PartitioningAlgorithm options: outsource searches to base templated
//!   on enum type
class PartitioningAlgorithm : public tk::Toggle< PartitioningAlgorithmType > {

  public:
    using ParamType = std::string;

    //! Valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::rcb
                                       , kw::rib
                                       , kw::hsfc
                                       , kw::mj
                                       , kw::phg
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PartitioningAlgorithm() :
      tk::Toggle< PartitioningAlgorithmType >(
        //! Group, i.e., options, name
        "Mesh partitioning algorithm",
        //! Enums -> names
        { { PartitioningAlgorithmType::RCB, kw::rcb::name() },
          { PartitioningAlgorithmType::RIB, kw::rib::name() },
          { PartitioningAlgorithmType::HSFC, kw::hsfc::name() },
          { PartitioningAlgorithmType::MJ, kw::mj::name() },
          { PartitioningAlgorithmType::PHG, kw::phg::name() } },
        //! keywords -> Enums
        { { kw::rcb::string(), PartitioningAlgorithmType::RCB },
          { kw::rib::string(), PartitioningAlgorithmType::RIB },
          { kw::hsfc::string(), PartitioningAlgorithmType::HSFC },
          { kw::mj::string(), PartitioningAlgorithmType::MJ },
          { kw::phg::string(), PartitioningAlgorithmType::PHG } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    const ParamType& param( PartitioningAlgorithmType m ) const {
      using tk::operator<<;
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for PartitioningAlgorithm \"")
              << m << "\"" );
      return it->second;
    }

    // Return true if partitioning algorithm is geometric
    //! \param[in] m Enum value of the option requested
    //! \return True if partitioning algorithm is geometric, false if it is not
    bool geometric( PartitioningAlgorithmType m ) const {
      if ( m == PartitioningAlgorithmType::RCB ||
           m == PartitioningAlgorithmType::RIB ||
           m == PartitioningAlgorithmType::HSFC ||
           m == PartitioningAlgorithmType::MJ )
        return true;
      else
       return false;
    }

  private:
    //! Enums -> Zoltan partitioning algorithm parameters
    std::map< PartitioningAlgorithmType, ParamType > method {
      { PartitioningAlgorithmType::RCB, "rcb" },
      { PartitioningAlgorithmType::RIB, "rib" },
      { PartitioningAlgorithmType::HSFC, "hsfc" },
      { PartitioningAlgorithmType::MJ, "multijagged" },
      { PartitioningAlgorithmType::PHG, "phg" }
    };
};

} // ctr::
} // tk::

#endif // InciterPartitioningAlgorithmOptions_h
