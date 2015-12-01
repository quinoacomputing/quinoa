//******************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.C
  \author    J. Bakosi
  \date      Tue 01 Dec 2015 09:32:49 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh graph
    partitioning.
*/
//******************************************************************************

#include <string>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <iosfwd>
#include <cstdlib>
#include <cstring>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Zoltan2_PartitioningSolution.hpp>

#include "Exception.h"
#include "UnsMesh.h"
#include "ExceptionMPI.h"
#include "ZoltanInterOp.h"
#include "DerivedData.h"
#include "Reorder.h"

namespace tk {
namespace zoltan {


//! GeometricMeshElemAdapter : Zoltan2::MeshAdapter
//! \details GeometricMeshElemAdapter specializes those virtual member functions
//!   of Zoltan2::MeshAdapter that are required for mesh-element-based
//!   geometric partitioning with Zoltan2
template< typename User >
class GeometricMeshElemAdapter : public Zoltan2::MeshAdapter< User > {

  private:
    using MeshEntityType = Zoltan2::MeshEntityType;
    using EntityTopologyType = Zoltan2::EntityTopologyType;

  public:
    using gno_t = typename Zoltan2::InputTraits< User >::gno_t;
    using scalar_t = typename Zoltan2::InputTraits< User >::scalar_t;
    using base_adapter_t = Zoltan2::MeshAdapter< User >;

    //! Constructor
    //! \param[in] nelem Number of elements in mesh graph on this rank
    //! \param[in] centroid Mesh element coordinates (centroids)
    //! \param[in] elemid Mesh element global IDs
    GeometricMeshElemAdapter(
      std::size_t nelem,
      const std::array< std::vector< tk::real >, 3 >& centroid,
      const std::vector< std::size_t >& elemid )
    : m_nelem( nelem ),
      m_topology( EntityTopologyType::TETRAHEDRON ),
      m_centroid( centroid ),
      m_elemid( elemid )
    {}

    //! Returns the number of mesh entities on this rank
    //! \return Number of mesh elements on this rank
    std::size_t getLocalNumOf( MeshEntityType etype ) const override
    { return m_nelem; }

    //! Provide a pointer to this rank's identifiers
    //! \param[inout] Ids Pointer to the list of global element Ids on this rank
    void getIDsViewOf( MeshEntityType etype, const gno_t*& Ids) const override
    { Ids = m_elemid.data(); }

    //! Provide a pointer to the entity topology types
    //! \param Types Pointer to the list of entity topology types on this rank
    void getTopologyViewOf( MeshEntityType etype,
                            const EntityTopologyType*& Types ) const override
    { Types = &m_topology; }

    //! Return dimensionality of the mesh
    //! \return Number of mesh dimension
    int getDimension() const override { return 3; }

    //! Provide a pointer to one dimension of mesh element coordinates
    //! \param[in] coords Pointer to a list of coordinate values for the
    //!   dimension
    //! \param[inout] stride Describes the layout of the coordinate values in
    //!   the coords list. If stride is one, then the ith coordinate value is
    //!   coords[i], but if stride is two, then the ith coordinate value is
    //!   coords[2*i]
    //! \param dim Value from 0 to one less than getEntityCoordinateDimension()
    //!   specifying which dimension is being provided in the coords list
    void getCoordinatesViewOf( MeshEntityType etype,
                               const scalar_t*& coords,
                               int &stride,
                               int dim ) const override
    {
      coords = m_centroid[ static_cast<std::size_t>(dim) ].data();
      stride = 1;
    }

  private:
    const std::size_t m_nelem;           //!< Number of elements on this rank
    const EntityTopologyType m_topology; //!< Mesh element topology types
    //! Mesh element coordinates (centroids)
    const std::array< std::vector< tk::real >, 3 >& m_centroid;
    //! Global mesh element ids
    const std::vector< std::size_t >& m_elemid;
};

std::vector< std::size_t >
geomPartMesh( tk::ctr::PartitioningAlgorithmType algorithm,
              const std::array< std::vector< tk::real >, 3 >& centroid,
              const std::vector< std::size_t >& elemid,
              std::size_t nelem,
              int npart )
//******************************************************************************
//  Partition mesh using Zoltan2 with a geometric partitioner, such as RCB, RIB
//! \param[in] algorithm Partitioning algorithm type
//! \param[in] centroid Mesh element coordinates
//! \param[in] elemid Global mesh element ids
//! \parampin] nelem Number of elements in mesh (on this MPI rank)
//! \param[in] npart Number of desired graph partitions
//! \return Array of chare ownership IDs mapping graph points to concurrent
//!   async chares
//! \details This function uses Zoltan to partition the mesh graph in parallel.
//!   It assumes that the mesh graph is distributed among all the MPI ranks.
//! \author J. Bakosi
//******************************************************************************
{
  // Set Zoltan parameters
  Teuchos::ParameterList params( "Zoltan parameters" );
  params.set( "algorithm", tk::ctr::PartitioningAlgorithm().param(algorithm) );
  params.set( "num_global_parts", std::to_string(npart) );
  params.set( "objects_to_partition", "mesh_elements" );

  // Define types for Zoltan2
  //  * 1st argument, 'scalar': the data type for element values, weights and
  //    coordinates
  //  * 2nd argument, 'gid' (global id): is the data type used by the
  //    application for global Ids. If the application's global Id data type is
  //    a Teuchos Ordinal, then gid and gno can the same. Otherwise, the
  //    application global Ids will be mapped to Teuchos Ordinals for use by
  //    Zoltan2 internally. (Teuchos Ordinals are those data types for which
  //    traits are defined in Teuchos_OrdinalTraits.hpp.)
  //  * 3rd argument, 'lno' (local number): the integral data type used by the
  //    application and by Zoltan2 for local indices and local counts
  //  * 4th argument 'gno' (global number): is the integral data type used by
  //    Zoltan2 to represent global indices and global counts
  using ZoltanTypes =
    Zoltan2::BasicUserTypes< tk::real, std::size_t, std::size_t >;

  // Create mesh adapter for Zoltan for mesh element partitioning
  using InciterZoltanAdapter = GeometricMeshElemAdapter< ZoltanTypes >;
  InciterZoltanAdapter adapter( nelem, centroid, elemid );

  // Create Zoltan2 partitioning problem using our mesh input adapter
  Zoltan2::PartitioningProblem< InciterZoltanAdapter >
    partitioner( &adapter, &params );

  // Perform partitioning using Zoltan
  partitioner.solve();

  // Copy over array of chare IDs corresponding to the ownership of elements
  // in our chunk of the mesh graph, i.e., the coloring or chare ids for the
  // mesh elements we operated on
  auto partlist = partitioner.getSolution().getPartListView();
  std::vector< std::size_t > chare( nelem );
  for (std::size_t p=0; p<nelem; ++p )
    chare[p] = static_cast< std::size_t >( partlist[p] );

  return chare;
}

} // zoltan::
} // tk::
