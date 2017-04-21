// *****************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh graph
    partitioning.
*/
// *****************************************************************************

#include "NoWarning/Zoltan2_MeshAdapter.h"
#include "NoWarning/Zoltan2_PartitioningProblem.h"
#include <Zoltan2_PartitioningSolution.hpp>

#include "ZoltanInterOp.h"

namespace tk {
namespace zoltan {

//! GeometricMeshElemAdapter : Zoltan2::MeshAdapter
//! \details GeometricMeshElemAdapter specializes those virtual member functions
//!   of Zoltan2::MeshAdapter that are required for mesh-element-based
//!   geometric partitioning with Zoltan2
template< typename ZoltanTypes >
class GeometricMeshElemAdapter : public Zoltan2::MeshAdapter< ZoltanTypes > {

  private:
    using MeshEntityType = Zoltan2::MeshEntityType;
    using EntityTopologyType = Zoltan2::EntityTopologyType;

  public:
    using gno_t = typename Zoltan2::InputTraits< ZoltanTypes >::gno_t;
    using scalar_t = typename Zoltan2::InputTraits< ZoltanTypes >::scalar_t;
    using base_adapter_t = Zoltan2::MeshAdapter< ZoltanTypes >;

    //! Constructor
    //! \param[in] nelem Number of elements in mesh graph on this rank
    //! \param[in] centroid Mesh element coordinates (centroids)
    //! \param[in] elemid Mesh element global IDs
    GeometricMeshElemAdapter(
      std::size_t nelem,
      const std::array< std::vector< tk::real >, 3 >& centroid,
      const std::vector< long >& elemid )
    : m_nelem( nelem ),
      m_topology( EntityTopologyType::TETRAHEDRON ),
      m_centroid( centroid ),
      m_elemid( elemid )
    {}

    //! Returns the number of mesh entities on this rank
    //! \return Number of mesh elements on this rank
    std::size_t getLocalNumOf( MeshEntityType ) const override
    { return m_nelem; }

    //! Provide a pointer to this rank's identifiers
    //! \param[in,out] Ids Pointer to the list of global element Ids on this
    //!   rank
    void getIDsViewOf( MeshEntityType, const gno_t*& Ids) const override
    { Ids = m_elemid.data(); }

    //! Provide a pointer to the entity topology types
    //! \param Types Pointer to the list of entity topology types on this rank
    void getTopologyViewOf( MeshEntityType,
                            const EntityTopologyType*& Types ) const override
    { Types = &m_topology; }

    //! Return dimensionality of the mesh
    //! \return Number of mesh dimension
    int getDimension() const override { return 3; }

    //! Provide a pointer to one dimension of mesh element coordinates
    //! \param[in] coords Pointer to a list of coordinate values for the
    //!   dimension
    //! \param[in,out] stride Describes the layout of the coordinate values in
    //!   the coords list. If stride is one, then the ith coordinate value is
    //!   coords[i], but if stride is two, then the ith coordinate value is
    //!   coords[2*i]
    //! \param dim Value from 0 to one less than getEntityCoordinateDimension()
    //!   specifying which dimension is being provided in the coords list
    void getCoordinatesViewOf( MeshEntityType,
                               const scalar_t*& coords,
                               int &stride,
                               int dim ) const override
    {
      coords = m_centroid[ static_cast<std::size_t>(dim) ].data();
      stride = 1;
    }

  private:
    //! Number of elements on this rank
    const std::size_t m_nelem;
    //! Mesh element topology types
    const EntityTopologyType m_topology;
    //! Mesh element coordinates (centroids)
    const std::array< std::vector< tk::real >, 3 >& m_centroid;
    //! Global mesh element ids
    const std::vector< long >& m_elemid;
};

std::vector< std::size_t >
geomPartMesh( tk::ctr::PartitioningAlgorithmType algorithm,
              const std::array< std::vector< tk::real >, 3 >& centroid,
              const std::vector< long >& elemid,
              std::size_t nelem,
              int npart )
// *****************************************************************************
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
// *****************************************************************************
{
  // Set Zoltan parameters
  Teuchos::ParameterList params( "Zoltan parameters" );
  params.set( "algorithm", tk::ctr::PartitioningAlgorithm().param(algorithm) );
  params.set( "num_global_parts", std::to_string(npart) );
  params.set( "objects_to_partition", "mesh_elements" );

  // Define types for Zoltan2
  //  * 1st argument, 'scalar': the data type for element values, weights and
  //    coordinates
  //  * 2nd argument, 'lno' (local number): the integral data type used by
  //    the application, i.e, quinoa, and by Zoltan2 for local indices and local
  //    counts
  //  * 3rd argument 'gno' (global number): is the integral data type used by
  //    the application, i.e., quinoa, and Zoltan2 to represent global
  //    identifiers and global counts
  // See also
  // tpl/src/trilinos/packages/zoltan2/src/input/Zoltan2_InputTraits.hpp
  using ZoltanTypes = Zoltan2::BasicUserTypes< tk::real, long, long >;

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
