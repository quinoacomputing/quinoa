// *****************************************************************************
/*!
  \file      src/Control/Options/PartitioningAlgorithm.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh partitioning algorithm options
  \details   Mesh partitioning algorithm options
*/
// *****************************************************************************
#ifndef InciterPartitioningAlgorithmOptions_h
#define InciterPartitioningAlgorithmOptions_h

#include <map>

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

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

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit PartitioningAlgorithm() :
      tk::Toggle< PartitioningAlgorithmType >(
        //! Group, i.e., options, name
        "Mesh partitioning algorithm",
        //! Enums -> names
        { { PartitioningAlgorithmType::RCB, "rcb" },
          { PartitioningAlgorithmType::RIB, "rib" },
          { PartitioningAlgorithmType::HSFC, "hsfc" },
          { PartitioningAlgorithmType::MJ, "mj" },
          { PartitioningAlgorithmType::PHG, "phg" } },
        //! keywords -> Enums
        { { "rcb", PartitioningAlgorithmType::RCB },
          { "rib", PartitioningAlgorithmType::RIB },
          { "hsfc", PartitioningAlgorithmType::HSFC },
          { "mj", PartitioningAlgorithmType::MJ },
          { "phg", PartitioningAlgorithmType::PHG } } ) {}

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
