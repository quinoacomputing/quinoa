// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/MeshVelocitySmoother.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh velocity smoother configuration options for inciter
  \details   Mesh velocity smoother configuration options for inciter.
*/
// *****************************************************************************
#ifndef MeshVelocitySmootherOptions_h
#define MeshVelocitySmootherOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Mesh velocity smoother configuration option types
enum class MeshVelocitySmootherType : uint8_t { NONE
                                              , LAPLACE
                                              , HELMHOLTZ
                                              };

//! \brief Pack/Unpack MeshVelocitySmootherType: forward overload to generic
//!    enum class packer
inline void operator|( PUP::er& p, MeshVelocitySmootherType& e )
{ PUP::pup( p, e ); }

//! \brief Mesh velocity options: outsource to base templated on enum type
class MeshVelocitySmoother : public tk::Toggle< MeshVelocitySmootherType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MeshVelocitySmoother() :
      tk::Toggle< MeshVelocitySmootherType >(
        //! Group, i.e., options, name
        "Mesh velocity smoother",
        //! Enums -> names (if defined, policy codes, if not, name)
        { { MeshVelocitySmootherType::NONE, "none" }
        , { MeshVelocitySmootherType::LAPLACE, "laplace" }
        , { MeshVelocitySmootherType::HELMHOLTZ, "helmholtz" }
        },
        //! keywords -> Enums
        { { "none", MeshVelocitySmootherType::NONE }
        , { "laplace", MeshVelocitySmootherType::LAPLACE }
        , { "helmholtz", MeshVelocitySmootherType::HELMHOLTZ }
        } )
    {}

};

} // ctr::
} // inciter::

#endif // MeshVelocitySmootherOptions_h
