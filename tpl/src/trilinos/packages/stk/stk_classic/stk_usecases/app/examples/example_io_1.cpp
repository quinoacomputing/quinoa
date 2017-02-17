/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>
#include <set>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace stk_examples {

typedef stk_classic::mesh::Field<double>                                ScalarField ;
typedef stk_classic::mesh::Field<int>                                   ScalarIntField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian>          CartesianField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::FullTensor>         FullTensorField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::SymmetricTensor>    SymmetricTensorField ;


CartesianField &
declare_vector_field_on_all_nodes(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  stk_classic::mesh::fem::FEMMetaData &fem = stk_classic::mesh::fem::FEMMetaData::get(meta_data);  
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), fem.node_rank() , meta_data.universal_part() , n1 );
}


CartesianField &
declare_vector_field_on_all_elements(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  
  stk_classic::mesh::fem::FEMMetaData &fem = stk_classic::mesh::fem::FEMMetaData::get(meta_data);  
  const stk_classic::mesh::EntityRank element_rank = fem.element_rank();

  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), element_rank , meta_data.universal_part() , n1 );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Example reader / writer

void example_io_1( stk_classic::ParallelMachine comm,
                   const std::string& in_filename,
                   const std::string& out_filename)
{
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  static const size_t spatial_dimension = 3;

  std::cout
 << "========================================================================\n"
 << " Use Case 1: Simple mesh I/O                                            \n"
 << "========================================================================\n";
 
  stk_classic::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk_classic::mesh::Part & universal        = meta_data.universal_part();
  stk_classic::mesh::put_field(meta_data.declare_field< CartesianField >( "coordinates" ) , meta_data.node_rank() , universal , spatial_dimension );

  //----------------------------------
  const std::string dbtype("exodusii");

  // Open, read, filter meta data from the input mesh file:
  // The coordinates field will be set to the correct dimension.

  stk_classic::io::MeshData mesh_data;
  stk_classic::io::create_input_mesh(dbtype, in_filename, comm, meta_data, mesh_data);

  //----------------------------------
  // Print the parts that were read from the file:

  const stk_classic::mesh::PartVector & io_parts = meta_data.get_parts();
  for (stk_classic::mesh::PartVector::const_iterator i = io_parts.begin();
           i != io_parts.end(); ++i) {
    stk_classic::mesh::Part & part = **i ;

    const CellTopologyData * const top = meta_data.get_cell_topology( part ).getCellTopologyData();

    const unsigned num_subsets = part.subsets().size();

    std::cout << "  Part: " << part.name();
    if ( part.primary_entity_rank() < meta_data.entity_rank_count() ) {
      std::cout << " "
                << meta_data.entity_rank_name( part.primary_entity_rank() );
    }
    if ( top ) { std::cout << " " << top->name ; }
    if ( num_subsets ) { std::cout << " with " << num_subsets << " subsets" ; }
    std::cout << std::endl ;
  }

  //----------------------------------
  // Now I can process what was read from the mesh file to
  // declare more parts and fields.

  meta_data.commit();

  //----------------------------------
  // Create mesh bulk data conforming to the mesh meta data
  // and distributed among the given parallel machine.

  stk_classic::mesh::BulkData bulk_data(meta_data.get_meta_data(meta_data), comm);

  // Read the model (topology, coordinates, attributes, etc)
  // from the mesh-file into the mesh bulk data.
  stk_classic::io::populate_bulk_data(bulk_data, mesh_data);

  //----------------------------------
  // Create a mesh writer that will simply write out what was read.
  // the parts, attributes, and transient arguments can be different
  // that what was read.

  stk_classic::io::create_output_mesh(out_filename, comm, bulk_data, mesh_data);

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace stk_examples


