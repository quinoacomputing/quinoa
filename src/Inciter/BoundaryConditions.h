// *****************************************************************************
/*!
  \file      src/Inciter/BoundaryConditions.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Data and functionality working on boundary conditions
  \details   Data and functionality working on boundary conditions.
*/
// *****************************************************************************
#ifndef BoundaryConditions_h
#define BoundaryConditions_h

#include <vector>
#include <map>
#include <unordered_map>

#include "Discretization.h"

#include "boundaryconditions.decl.h"

namespace inciter {

//! Data and functionality working on boundary conditions
class BoundaryConditions : public CBase_BoundaryConditions {

  public:
    //! Constructor
    explicit BoundaryConditions(
      const std::map< int, std::vector< std::size_t > >& sidenodes );

    //! Create map that assigns the local mesh node IDs mapped to side set ids
    std::map< int, std::vector< std::size_t > >
    sideNodes( const std::unordered_map< std::size_t, std::size_t >& filenodes,
               const std::unordered_map< std::size_t, std::size_t >& lid );

  private:
    //! Map associating file-node IDs to side set IDs for all side sets in file
    //! \details This map stores mesh node IDs as exist in the mesh file
    //!   associated to side set IDs for all side sets read from mesh file
    //!   independent of what the user sets boundary conditions on
    std::map< int, std::vector< std::size_t > > m_sideFileNodes;
};

} // inciter::

#endif // BoundaryConditions_h
