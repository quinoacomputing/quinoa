#include "Mesquite_MeshTransform.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_CLArgs.hpp"
#include "Mesquite_MsqError.hpp"

#include <iostream>

using namespace Mesquite;

const char ROTATE_FLAG = 'r';
const char SCALE_FLAG = 's';
const char TRANSLATE_FLAG = 't';

class RotateArg : public CLArgs::DoubleListArgI
{
  private: MeshTransform* mTransform;
  public:  RotateArg( MeshTransform* xform ) : mTransform(xform) {}
      bool value( const std::vector<double>& vals );
};
class ScaleArg : public CLArgs::DoubleListArgI
{
  private: MeshTransform* mTransform;
  public:  ScaleArg( MeshTransform* xform ) : mTransform(xform) {}
      bool value( const std::vector<double>& vals );
};
class TranslateArg : public CLArgs::DoubleListArgI
{
  private: MeshTransform* mTransform;
  public:  TranslateArg( MeshTransform* xform ) : mTransform(xform) {}
      bool value( const std::vector<double>& vals );
};
bool RotateArg::value( const std::vector<double>& vals )
{
  mTransform->add_rotation( Vector3D(vals[1],vals[2],vals[3]), vals[0] * M_PI / 180 );
  return true;
}
bool ScaleArg::value( const std::vector<double>& vals )
{
  for (unsigned i = 0; i < vals.size(); ++i)
    if (vals[i] <= 0)
      return false;
  if (vals.size() == 1)
    mTransform->add_scale( vals[0] );
  else
    mTransform->add_scale( Vector3D(vals[0], vals[1], vals[2]) );
  return true;
}
bool TranslateArg::value( const std::vector<double>& vals )
{
  mTransform->add_translation( Vector3D(vals[0], vals[1], vals[2]) );
  return true;
}

int main( int argc, char* argv[] )
{
  MeshTransform xform;
  RotateArg rotate_arg( &xform );
  ScaleArg scale_arg( &xform );
  TranslateArg translate_arg( &xform );
  CLArgs::ToggleArg freeonly, skin;

  CLArgs args( "vtkxform", "Transform a mesh",
               "Apply one or more transformations to vertex coordinates "
               "in a mesh read from a VTK file." );
  
  const char* ROTATE_VALS[] = { "a", "i", "j", "k" };
  args.double_list_flag( ROTATE_FLAG, "Specify a rotation as an angle in degrees counter-clockwise about a vector", &rotate_arg );
  args.limit_list_flag( ROTATE_FLAG, 4, ROTATE_VALS );
  
  const char* SCALE_VALS[] = { "s", "sx", "sy", "sz" };
  args.double_list_flag( SCALE_FLAG, "Specify factor(s) by which to scale mesh about origin", &scale_arg );
  args.limit_list_flag( SCALE_FLAG, 1, SCALE_VALS );
  args.limit_list_flag( SCALE_FLAG, 3, SCALE_VALS + 1 );
  
  const char* TRANSLATE_VALS[] = { "dx", "dy", "dz" };
  args.double_list_flag( TRANSLATE_FLAG, "Specify translation of vertex coordinates.", &translate_arg );
  args.limit_list_flag( TRANSLATE_FLAG, 3, TRANSLATE_VALS );

  args.toggle_flag( 'f', "Do not move fixed vertices.", &freeonly );
  args.toggle_flag( 'k', "Mark boundary vertices as fixed", &skin );

  args.add_required_arg( "input_file" );
  args.add_required_arg( "output_file" );
  
  
  std::vector<std::string> files;
  if (!args.parse_options( argc, argv, files, std::cerr )) {
    args.print_usage( std::cerr );
    exit(1);
  }
  std::string input_file = files[0];
  std::string output_file = files[1];

  MeshImpl mesh;
  MsqError err;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << err << std::endl 
              << "Failed to read file: " << input_file << std::endl;
    return 1;
  }
  
  if (skin.value()) {
    mesh.mark_skin_fixed( err, false );
    if (err) {
      std::cerr << err << std::endl
                << "Failed to skin mesh from file: " << input_file << std::endl;
      return 1;
    }
  }
  
  xform.skip_fixed_vertices( freeonly.value() );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, 0);
  xform.loop_over_mesh( &mesh_and_domain, 0, err );
  if (err) {
    std::cerr << err << std::endl ;
    return 2;
  }
  
  mesh.write_vtk( output_file.c_str(), err );
  if (err) {
    std::cerr << err << std::endl 
                    << "Failed to write file: " << output_file << std::endl;
    return 1;
  }
  
  return 0;
}

