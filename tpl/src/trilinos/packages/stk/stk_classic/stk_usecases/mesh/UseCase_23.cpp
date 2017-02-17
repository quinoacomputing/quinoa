/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <stdexcept>
#include <sstream>
#include <vector>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Property.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <common/gnu_malloc_hooks.hpp>
#include <mesh/UseCase_14.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 23: internal force computation with
// a demonstration of part properties.
// The function 'use_case_23_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk_classic {
namespace app {

enum { SpatialDim = 3 };

namespace {

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::fem::FEMMetaData& meta_data)
{
  return meta_data.element_rank();
}

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::Part& part)
{
  return get_element_rank(stk_classic::mesh::fem::FEMMetaData::get(part));
}

CartesianField &
declare_vector_field_on_all_nodes(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), stk_classic::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() , n1 );
}


CartesianField &
declare_vector_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


ScalarField &
declare_scalar_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<ScalarField>(s), get_element_rank(meta_data), meta_data.universal_part() );
}


SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return put_field( meta_data.declare_field<SymmetricTensorField>(s), get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> &
put_field_on_elements( stk_classic::mesh::Field<Type,T1> & f , stk_classic::mesh::Part & p , unsigned n1 )
{
  stk_classic::mesh::put_field( f , get_element_rank(p) , p , n1 );
  return f ;
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> & put_field_on_all_elements( stk_classic::mesh::Field<Type,T1> & f , unsigned n1 )
{
  put_field_on_elements( f , stk_classic::mesh::fem::FEMMetaData::get(f).universal_part() , n1 );
  return f ;
}

} // namespace <unnamed>

void use_case_23_driver(
  MPI_Comm comm ,
  const std::string & file_name ,
  const unsigned box_size[] ,
  const unsigned box_sides[][2] ,
  const unsigned num_trials );

//--------------------------------------------------------------------
//
// main driver for use-case 23: element internal force with Property usage.
//

void use_case_23_driver( MPI_Comm comm , bool performance_test )
{
  int num_procs = stk_classic::parallel_machine_size( comm );

  if ( ! stk_classic::parallel_machine_rank( comm ) ) {
    std::cout << " stk_mesh Use Case #23 - element internal force with part properties, begin" << std::endl ;
  }

  if ( ! performance_test ) {
    // Quick & small test for correctness:
    const unsigned a_box_size[3] = { 10 , 10 , 10*num_procs };
    const unsigned a_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    int num_trials = 1;
    use_case_23_driver( comm , "use_case_23.exo" , a_box_size , a_box_sides , num_trials );
  }
  else {
    int num_trials = 582 ; // 582 ;

    const unsigned a_box_size[3] = { 100 , 100 , 100 };
    const unsigned a_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    std::cout << "Running 100x100x100 case, num_trials = " << num_trials << std::endl;

    use_case_23_driver( comm , "use_case_23a.exo" , a_box_size , a_box_sides , num_trials );

    const unsigned b_box_size[3] = { 100 , 125 , 10*num_procs };
    const unsigned b_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    std::cout << "Running 100x125x10*nprocs case, num_trials = " << num_trials << std::endl;

    use_case_23_driver( comm , "use_case_23b.exo" , b_box_size , b_box_sides , num_trials );
  }
}

//--------------------------------------------------------------------

inline
void zero_data( double * beg , double * const end )
{ while ( beg < end ) { *beg++ = 0 ; } }

template< class FieldType >
void zero_field_data( stk_classic::mesh::BulkData & mesh , unsigned type , const FieldType & field )
{
  typedef stk_classic::mesh::BucketArray< FieldType > array_type ;

  const std::vector<stk_classic::mesh::Bucket*> & ks = mesh.buckets( type );

  for ( std::vector<stk_classic::mesh::Bucket*>::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    stk_classic::mesh::Bucket & bucket = **ik ;
    array_type data( field , bucket );
    zero_data( data.contiguous_data() , data.contiguous_data() + data.size() );
  }
}

//--------------------------------------------------------------------


void use_case_23_driver(
  MPI_Comm comm ,
  const std::string & /*file_name */,
  const unsigned box_size[] ,
  const unsigned box_sides[][2] ,
  const unsigned num_trials )
{
  //----------------------------------

  APSHex8ug hex_element;

  lame::matParams materialParameters;

  // Timing:
  //   [0] = stk_classic::mesh::MetaData creation
  //   [1] = stk_classic::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force
  //   [4] = Parallel swap-add of internal force

  double time_min[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double time_max[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double wtime = 0 ;

  //--------------------------------------------------------------------

  reset_malloc_stats();

  if ( 0 == stk_classic::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh performance use case #23" << std::endl
              << "  Number Processes = " << stk_classic::parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  {
    wtime = stk_classic::wall_time();

    //------------------------------------------------------------------
    // Declare the mesh meta data and bulk data.

    stk_classic::mesh::fem::FEMMetaData mesh_meta_data( SpatialDim );
    const stk_classic::mesh::EntityRank element_rank = mesh_meta_data.element_rank();
    stk_classic::mesh::BulkData mesh_bulk_data( mesh_meta_data.get_meta_data(mesh_meta_data), MPI_COMM_WORLD, 1000 );

    //--------------------------------
    // Element-block declarations typically occur when reading the
    // mesh-file meta-data, and thus won't usually appear in application code.
    // Declaring the element blocks and associating an element traits
    // with each element block.

    stk_classic::mesh::Part & block_hex = mesh_meta_data.declare_part("block_1", element_rank);
    stk_classic::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
    stk_classic::mesh::fem::set_cell_topology( block_hex, hex_top );

    double dt  = 1.0e-02;
    double YM  = 1e7; // Young's Modulus
    double PR  = 0.33; // Poisson Ratio
    double element_mass_dummy_value = 1.0;
    std::string material_model_name("ELASTIC");

    double current_stable_time_step = dt;


    //------------------------
    //store dt on the mesh_meta_data as a Part property on the part 'block_hex'

    stk_classic::mesh::Property<double>& delta_t = mesh_meta_data.get_meta_data(mesh_meta_data).declare_property<double>("dt");

    mesh_meta_data.get_meta_data(mesh_meta_data).put_property( delta_t, block_hex );

    double* delta_t_ptr = stk_classic::mesh::property_data(delta_t, block_hex);

    *delta_t_ptr = dt;


    //------------------------

    //----------------------------------
    // Compute other material constants.
    double BM = YM/(3.0-6.0*PR); // Bulk Modulus
    double SM = YM/(2.0+2.0*PR); // Shear Modulus
    double TM = 2.0 * SM;
    double LAMBDA = YM*PR/((1.0+PR)*(1.0-2.0*PR));
    double DM = TM + LAMBDA; // Dilatational Modulus

    lame::MatProps materialProperties;

    std::vector<double> youngs_modulus; youngs_modulus.push_back(YM);
    std::vector<double> poissons_ratio; poissons_ratio.push_back(PR);
    std::vector<double> bulk_modulus;     bulk_modulus.push_back(BM);
    std::vector<double> shear_modulus;   shear_modulus.push_back(SM);
    std::vector<double> two_mu;                 two_mu.push_back(TM);
    std::vector<double> lambda;                 lambda.push_back(LAMBDA);
    std::vector<double> dilatational_modulus; dilatational_modulus.push_back(DM);

    materialProperties["YOUNGS_MODULUS"] = youngs_modulus;
    materialProperties["POISSONS_RATIO"] = poissons_ratio;
    materialProperties["BULK_MODULUS"]   = bulk_modulus;
    materialProperties["SHEAR_MODULUS"]  = shear_modulus;
    materialProperties["TWO_MU"]         = two_mu;
    materialProperties["LAMBDA"]         = lambda;
    materialProperties["DILATATIONAL_MODULUS"] = dilatational_modulus;

    lame::Material * matmodel = lame::Elastic::createMaterial( materialProperties );


    //------------------------
    //store materialProperties on the mesh_meta_data as a Part property on the part 'block_hex'

    stk_classic::mesh::Property<lame::MatProps>& mprops = mesh_meta_data.get_meta_data(mesh_meta_data).declare_property<lame::MatProps>("materialProperties");
    mesh_meta_data.get_meta_data(mesh_meta_data).put_property(mprops, block_hex);

    lame::MatProps* mprops_ptr = stk_classic::mesh::property_data(mprops, block_hex);

    *mprops_ptr = materialProperties;


    //--------------------------------

    // Nodal vector fields
    CartesianField &model_coordinates     = declare_vector_field_on_all_nodes( mesh_meta_data , "model_coordinates" , SpatialDim );
    CartesianField &coordinates_field     = declare_vector_field_on_all_nodes( mesh_meta_data , "coordinates" , SpatialDim );
    CartesianField &velocity_field        = declare_vector_field_on_all_nodes( mesh_meta_data , "velocity" ,    SpatialDim );
    CartesianField &fint_field           = declare_vector_field_on_all_nodes( mesh_meta_data , "force_internal" , SpatialDim );

    // Element vector fields:
    CartesianField &Vorticity             = declare_vector_field_on_all_elements( mesh_meta_data , "Vorticity" , SpatialDim );

    // Element scalar fields:
    ScalarField &Shear_Modulus         = declare_scalar_field_on_all_elements( mesh_meta_data , "shear_modulus" );
    ScalarField &Dilatational_Modulus  = declare_scalar_field_on_all_elements( mesh_meta_data , "dilatational_modulus");
    ScalarField &Material_eff_twomu    = declare_scalar_field_on_all_elements( mesh_meta_data , "material_effictive_two_mu");
    ScalarField &Material_eff_bulk_mod = declare_scalar_field_on_all_elements( mesh_meta_data , "material_effective_bulk_moduli");
    ScalarField &Midstep_volume        = declare_scalar_field_on_all_elements( mesh_meta_data , "mid_step_volume");
    ScalarField &Element_time_step          = declare_scalar_field_on_all_elements( mesh_meta_data , "element_time_step");
    ScalarField &Element_mass          = declare_scalar_field_on_all_elements( mesh_meta_data , "element_mass");
    ScalarField &Hourglass_energy      = declare_scalar_field_on_all_elements( mesh_meta_data , "hourglass_energy");
    ScalarField &Internal_energy       = declare_scalar_field_on_all_elements( mesh_meta_data , "internal_energy");

    // Element symmetric tensor fields:
    SymmetricTensorField &Stretch       = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "stretch" ,    6 );
    SymmetricTensorField &StrainRate    = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "StrainRate" , 6 );
    SymmetricTensorField &RotatedStress = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "RotatedStress" , 6 );

    //--------------------------------
    // We don't like specifying the template type TWICE! very error prone.
    //--------------------------------
    // The multi-state fields don't have the 'declare and put' convenience functions (yet)
    //
    //  For clarity declare a integer to used for the number of states arguments.
    //
    const unsigned two_states = 2 ;

    // Element two state symmetric tensor field, on all elements (as two function calls):
    SymmetricTensorField &StressNew = mesh_meta_data.declare_field< SymmetricTensorField >( "Stress" , two_states );
    put_field_on_all_elements( StressNew , 6 );

    // Element two state full tensor field on all elements (as nested function calls):
    FullTensorField &RotationNew =
      put_field_on_all_elements( mesh_meta_data.declare_field< FullTensorField >("Rotation", two_states ) , 9 );

    //--------------------------------
    // Hourglass fields, these don't have the convenience functions for 'declare and put'

    HourglassArrayField &HourglassResistanceNew =
      put_field_on_all_elements( mesh_meta_data.declare_field< HourglassArrayField >("HourglassResistance" , two_states ) , 12 );

    HourglassOpField & MidHourglassOp =
      put_field_on_all_elements( mesh_meta_data.declare_field< HourglassOpField >("mid_hourglass_operator") , 32 );

    //--------------------------------
    // Declare aggressive "gather" fields which are an array of
    // pointers to the element's nodes' coordinate, velocity, and
    // internal force field data.
    //
    // The declarations specify element fields of the following form:
    //
    //     double * coord_gather   [ nodes_per_element ]
    //     double * velocity_gather[ nodes_per_element ]
    //     double * fint_gather    [ nodes_per_element ]
    //
    // where
    //
    //     coord_gather[i]    == field_data( coordinates_field , element_node[i] )
    //     velocity_gather[i] == field_data( velocity ,          element_node[i] )
    //     fint_gather[i]     == field_data( fint ,              element_node[i] )
    //
    // The number of nodes per element could vary, so the field is put on each element block
    // with a size of the number of nodes per element in that element block.

    ElementNodePointerField & coord_gather    = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "coord_gather" , coordinates_field );
    ElementNodePointerField & velocity_gather = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "velocity_gather" , velocity_field );
    ElementNodePointerField & fint_gather     = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "fint_gather" , fint_field );

    put_field_on_elements( coord_gather ,    block_hex , shards::Hexahedron<> ::node_count );
    put_field_on_elements( velocity_gather , block_hex , shards::Hexahedron<> ::node_count );
    put_field_on_elements( fint_gather ,     block_hex , shards::Hexahedron<> ::node_count );

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    mesh_meta_data.commit();

    //------------------------------------------------------------------

    time_max[0] = stk_classic::wall_dtime( wtime );

    //------------------------------------------------------------------


    // In a typical app, the mesh would be read from file at this point.
    // But in this use-case, we generate the mesh and initialize
    // field data to use-case defined values.

    stk_classic::app::use_case_14_generate_mesh(
      mesh_bulk_data ,
      box_size ,
      model_coordinates ,
      coord_gather ,
      block_hex ,
      box_sides );

    stk_classic::app::use_case_14_initialize_data(
      mesh_bulk_data ,
      model_coordinates ,
      coordinates_field ,
      velocity_field );

    time_max[1] = stk_classic::wall_dtime( wtime );

  //------------------------------------------------------------------
  // Ready to run the algorithms:
  //------------------------------------------------------------------


  stk_classic::mesh::Selector select_owned( stk_classic::mesh::MetaData::get(mesh_bulk_data).locally_owned_part() );

  const std::vector< stk_classic::mesh::Bucket * > & element_buckets =
    mesh_bulk_data.buckets(element_rank);

  unsigned maximum_bucket_size = 0 ;
  for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
    const stk_classic::mesh::Bucket & bucket = **k ;
    if ( maximum_bucket_size < bucket.size() ) { maximum_bucket_size = bucket.size(); }
  }

  // Need both the the old and new states of these two-state fields:

  HourglassArrayField  & HourglassResistanceOld = HourglassResistanceNew.field_of_state( stk_classic::mesh::StateOld );
  SymmetricTensorField & StressOld              = StressNew.field_of_state( stk_classic::mesh::StateOld );
  FullTensorField      & RotationOld            = RotationNew.field_of_state( stk_classic::mesh::StateOld );

  for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
    const stk_classic::mesh::Bucket & bucket = **k ;

    const int num_elements = bucket.size();

    double * stretch        = stk_classic::mesh::field_data( Stretch,    bucket.begin() );
    double * strain_rate    = stk_classic::mesh::field_data( StrainRate, bucket.begin() );
    double * stress_new     = stk_classic::mesh::field_data( StressNew,  bucket.begin() );
    double * stress_old     = stk_classic::mesh::field_data( StressOld , bucket.begin() );
    double * rotated_stress = stk_classic::mesh::field_data( RotatedStress , bucket.begin() );
    double * rotation_old   = stk_classic::mesh::field_data( RotationOld , bucket.begin() );
    double * rotation_new   = stk_classic::mesh::field_data( RotationNew , bucket.begin() );
    double * mass           = stk_classic::mesh::field_data( Element_mass , bucket.begin() );
    double * hg_old         = stk_classic::mesh::field_data( HourglassResistanceOld , bucket.begin() );
    double * hg_new         = stk_classic::mesh::field_data( HourglassResistanceNew , bucket.begin() );
    double *vorticity_ptr   = stk_classic::mesh::field_data( Vorticity, bucket.begin());
    double *mid_hg_op_ptr   = stk_classic::mesh::field_data( MidHourglassOp, bucket.begin());

    for ( int i = 0 ; i < num_elements ; ++i ) {

      stk_classic::Copy<32>( mid_hg_op_ptr , 0.0 ); mid_hg_op_ptr += 32 ;
      stk_classic::Copy< 3>( vorticity_ptr , 0.0 ); vorticity_ptr += 3 ;
      stk_classic::Copy< 6>( rotated_stress, 0.0 ); rotated_stress += 6 ;
      stk_classic::Copy< 6>( strain_rate,    0.0 ); strain_rate += 6 ;
      stk_classic::Copy<12>( hg_old,         0.0 ); hg_old += 12 ;
      stk_classic::Copy<12>( hg_new,         0.0 ); hg_new += 12 ;
      stk_classic::Copy< 6>( stress_new,     0.0 ); stress_new += 6 ;
      stk_classic::Copy< 6>( stress_old,     0.0 ); stress_old += 6 ;

      mass[i] = element_mass_dummy_value;

      // initialize stretch to identity.
      stretch[0] = 1.0;
      stretch[1] = 1.0;
      stretch[2] = 1.0;
      stretch[3] = 0.0;
      stretch[4] = 0.0;
      stretch[5] = 0.0;
      stretch+=6;

      // initialize rotation to identity.
      rotation_old[0] = 1.0;
      rotation_old[1] = 1.0;
      rotation_old[2] = 1.0;
      rotation_old[3] = 0.0;
      rotation_old[4] = 0.0;
      rotation_old[5] = 0.0;
      rotation_old[6] = 0.0;
      rotation_old[7] = 0.0;
      rotation_old[8] = 0.0;
      rotation_old +=9;

      // initialize rotation to identity.
      rotation_new[0] = 1.0;
      rotation_new[1] = 1.0;
      rotation_new[2] = 1.0;
      rotation_new[3] = 0.0;
      rotation_new[4] = 0.0;
      rotation_new[5] = 0.0;
      rotation_new[6] = 0.0;
      rotation_new[7] = 0.0;
      rotation_new[8] = 0.0;
      rotation_new += 9;
    }
  }
  std::cout << "KHP: Done element field init\n";

  for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
    const stk_classic::mesh::Bucket & bucket = **k ;

    const unsigned num_elements = bucket.size();

    double * const sm = stk_classic::mesh::field_data( Shear_Modulus , bucket.begin() );
    double * const dm = stk_classic::mesh::field_data( Dilatational_Modulus , bucket.begin() );
    double * const twomu = stk_classic::mesh::field_data( Material_eff_twomu , bucket.begin() );
    double * const bulk = stk_classic::mesh::field_data( Material_eff_bulk_mod , bucket.begin() );
    double * const mv = stk_classic::mesh::field_data( Midstep_volume , bucket.begin() );
    double * const hg_energy = stk_classic::mesh::field_data( Hourglass_energy , bucket.begin() );
    double * const internal_energy = stk_classic::mesh::field_data( Internal_energy , bucket.begin() );

    for ( unsigned i = 0 ; i < num_elements ; ++i ) {
      sm[i]              = SM;
      dm[i]              = DM;
      twomu[i]           = TM;
      bulk[i]            = BM;
      mv[i]              = 0.0;
      hg_energy[i]       = 0.0;
      internal_energy[i] = 0.0;
    }
  }

  //------------------------------------------------------------------
    time_max[2] = stk_classic::wall_dtime( wtime );
  //------------------------------------------------------------------

  // Scratch space

  enum { num_dof = 24 };
  const int num_dof_max_bucket = num_dof * maximum_bucket_size ;
  std::vector< double > vel(                 num_dof_max_bucket );
  std::vector< double > element_coordinates( num_dof_max_bucket );
  std::vector< double > force_new(           num_dof_max_bucket );

  wtime = stk_classic::wall_time();

  for(unsigned n=0; n<num_trials; ++n) {
    //
    // Call Internal Force!!!
    //

    // Need to zero out the old accumulated internal force so that its does not pollute the new accumulation.
    zero_field_data( mesh_bulk_data , stk_classic::mesh::fem::FEMMetaData::NODE_RANK , fint_field );

    for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
          k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
      const stk_classic::mesh::Bucket & bucket = **k ;

      const int num_elements = bucket.size();
      double *mid_hg_op_ptr = stk_classic::mesh::field_data( MidHourglassOp, bucket.begin());
      double *material_eff_twomu_ptr = stk_classic::mesh::field_data( Material_eff_twomu, bucket.begin());
      double *material_eff_bulk_modulus_ptr = stk_classic::mesh::field_data( Material_eff_bulk_mod, bucket.begin());
      double *mid_step_volume_ptr = stk_classic::mesh::field_data( Midstep_volume, bucket.begin());
      double *element_time_step_ptr = stk_classic::mesh::field_data( Element_time_step, bucket.begin());
      double *element_mass_ptr = stk_classic::mesh::field_data( Element_mass, bucket.begin());
      double *hg_energy_ptr = stk_classic::mesh::field_data( Hourglass_energy, bucket.begin());
      double *internal_energy_ptr = stk_classic::mesh::field_data (Internal_energy, bucket.begin());
      double *shear_modulus_ptr =  stk_classic::mesh::field_data( Shear_Modulus, bucket.begin());
      double *dilatational_modulus_ptr =  stk_classic::mesh::field_data( Dilatational_Modulus, bucket.begin() );

      double *rotation_old_ptr = stk_classic::mesh::field_data( RotationOld, bucket.begin());
      double *rotation_new_ptr = stk_classic::mesh::field_data( RotationNew, bucket.begin());

      double *stretch_ptr      = stk_classic::mesh::field_data( Stretch, bucket.begin());
      double *strain_rate_ptr  = stk_classic::mesh::field_data( StrainRate, bucket.begin());
      double *stress_old_ptr   = stk_classic::mesh::field_data( StressOld, bucket.begin());
      double *stress_new_ptr   = stk_classic::mesh::field_data( StressNew, bucket.begin());
      double *rotated_stress_ptr =  stk_classic::mesh::field_data( RotatedStress, bucket.begin());

      double *hg_resistance_old_ptr =  stk_classic::mesh::field_data( HourglassResistanceOld, bucket.begin());
      double *hg_resistance_new_ptr =  stk_classic::mesh::field_data( HourglassResistanceNew, bucket.begin());

      double *vorticity_ptr    =  stk_classic::mesh::field_data( Vorticity, bucket.begin());

      const double **coord    =
        (const double **) stk_classic::mesh::field_data( coord_gather,    bucket.begin());
      const double **velocity =
        (const double **) stk_classic::mesh::field_data( velocity_gather, bucket.begin());
      double **fint        = stk_classic::mesh::field_data( fint_gather, bucket.begin());

      //--------------------------------
      { // Gather nodal data into contiguous arrays.

        // element-node pointer fields
        const double ** field_coord  = coord ;
        const double ** vnodes       = velocity ;

        // scratch space for gather:
        double * elem_coord = & element_coordinates[0] ;
        double * elem_vel   = & vel[0] ;

        const int num_elem_nodes = num_elements * 8 ;

        for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
          const double * const f_coord = *field_coord ;
          const double * const f_vel   = *vnodes ;

          elem_coord[0] = f_coord[0];
          elem_coord[1] = f_coord[1];
          elem_coord[2] = f_coord[2];

          elem_vel[0] = f_vel[0];
          elem_vel[1] = f_vel[1];
          elem_vel[2] = f_vel[2];

          ++field_coord ;
          ++vnodes ;
          elem_coord += 3 ;
          elem_vel += 3 ;
        }
      }

      //--------------------------------

      double* dt_ptr = stk_classic::mesh::property_data(delta_t, block_hex);
      materialParameters.dt          = *dt_ptr;
      materialParameters.nelements   = num_elements ;
      materialParameters.strain_rate = strain_rate_ptr;
      materialParameters.stress_old  = stress_old_ptr;
      materialParameters.stress_new  = stress_new_ptr;

      lame::MatProps* materialProperties_ptr = stk_classic::mesh::property_data(mprops, block_hex);

      hex_element.internalForce( num_elements ,
        *dt_ptr,
        current_stable_time_step,
        element_time_step_ptr,
        *matmodel,
        materialParameters,
        *materialProperties_ptr,
        & element_coordinates[0],
        & vel[0],
        rotation_old_ptr,
        rotation_new_ptr,
        mid_step_volume_ptr,
        vorticity_ptr,
        stretch_ptr,
        strain_rate_ptr,
        mid_hg_op_ptr,
        stress_old_ptr,
        stress_new_ptr,
        rotated_stress_ptr,
        material_eff_bulk_modulus_ptr,
        material_eff_twomu_ptr,
        shear_modulus_ptr,
        dilatational_modulus_ptr,
        element_mass_ptr,
        & force_new[0],
        hg_energy_ptr,
        internal_energy_ptr,
        hg_resistance_old_ptr,
        hg_resistance_new_ptr
      );

      { // Scatter internal force values.
        double ** f = fint ;
        const double * f_new = & force_new[0] ;
        const int num_elem_nodes = num_elements * 8 ;
        for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
          double * const node_f = *f ;
          node_f[0] += f_new[0];
          node_f[1] += f_new[1];
          node_f[2] += f_new[2];
          ++f ;
          f_new += 3 ;
        }
      }
    }

    time_max[3] += stk_classic::wall_dtime( wtime );

    stk_classic::mesh::parallel_reduce( mesh_bulk_data , stk_classic::mesh::sum(fint_field) );

    time_max[4] += stk_classic::wall_dtime( wtime );

  }//end for(..num_trials...

  //------------------------------------------------------------------

  delete matmodel;
#ifdef USE_GNU_MALLOC_HOOKS
  if (parallel_machine_rank(comm) == 0) {
    double net_alloc = alloc_MB() - freed_MB();
    std::cout << "Mesh creation:" << "\n   Total allocated: "
       << alloc_MB()<<"MB in "<<alloc_blks() << " blocks."
       << "\n   Total freed: " << freed_MB() << "MB in "
       << freed_blks() << " blocks."
       << "\n   Net allocated: "<<net_alloc << "MB."<<std::endl;
  }
#endif

    //------------------------------------------------------------------
  }

  time_max[8] = stk_classic::wall_dtime( wtime );

  time_min[0] = time_max[0] ;
  time_min[1] = time_max[1] ;
  time_min[2] = time_max[2] ;
  time_min[3] = time_max[3] ;
  time_min[4] = time_max[4] ;
  time_min[5] = time_max[5] ;
  time_min[6] = time_max[6] ;
  time_min[7] = time_max[7] ;
  time_min[8] = time_max[8] ;

  stk_classic::all_reduce( comm , stk_classic::ReduceMax<9>( time_max ) & stk_classic::ReduceMin<9>( time_min ) );

  time_max[3] /= num_trials ;
  time_max[4] /= num_trials ;
  time_max[5] /= num_trials ;
  time_max[6] /= num_trials ;

  time_min[3] /= num_trials ;
  time_min[4] /= num_trials ;
  time_min[5] /= num_trials ;
  time_min[6] /= num_trials ;

  //   [0] = stk_classic::mesh::MetaData creation
  //   [1] = stk_classic::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force

  if ( ! stk_classic::parallel_machine_rank( comm ) ) {
    std::cout
      << "stk_mesh performance use case results:" << std::endl
      << "  Number trials       = " << num_trials << std::endl
      << "  Box size            = { " << box_size[0] << " , "
                                      << box_size[1] << " , "
                                      << box_size[2] << " }" << std::endl
      << "  Box sides           = { { "
                                    << box_sides[0][0] << " , "
                                    << box_sides[0][1] << " } , { "
                                    << box_sides[1][0] << " , "
                                    << box_sides[1][1] << " } , { "
                                    << box_sides[2][0] << " , "
                                    << box_sides[2][1] << " } }"
                                    << std::endl
      << "  Meta-data setup     = " << time_min[0] << " : "
                                    << time_max[0] << " sec, min : max"
                                    << std::endl
      << "  Bulk-data generation= " << time_min[1] << " : "
                                    << time_max[1] << " sec, min : max"
                                    << std::endl
      << "  Initialization      = " << time_min[2] << " : "
                                    << time_max[2] << " sec, min : max"
                                    << std::endl
      << "  Internal force      = " << time_min[3] << " : "
                                    << time_max[3] << " sec, min : max"
                                    << std::endl
      << "  Internal force (total) = " << time_min[3]*num_trials
                                    << std::endl
      << "  Swap-add            = " << time_min[4] << " : "
                                    << time_max[4] << " sec, min : max"
                                    << std::endl
      << "  Mesh destruction    = " << time_min[8] << " : "
                                    << time_max[8] << " sec, min : max"
                                    << std::endl
      << std::endl ;
  }
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace app
} // namespace stk_classic

