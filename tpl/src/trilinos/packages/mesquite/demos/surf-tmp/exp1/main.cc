#include "mesquite_version.h"

#include "TMPQualityMetric.hpp"
#if MSQ_VERSION_MINOR == 99
#  include "IdealShapeTarget.hpp"
#else
#  include "IdealTargetCalculator.hpp"
   typedef Mesquite::IdealTargetCalculator IdealShapeTarget;
#endif

#include "Target2DShapeBarrier.hpp"
#include "Target2DShape.hpp"
typedef Mesquite::Target2DShape TMetric;
#define TMetric_NAME "Shape"

#include "QualityAssessor.hpp"

#include "PlanarDomain.hpp"
#include "SphericalDomain.hpp"
#include "CylinderDomain.hpp"
#include "ConicDomain.hpp"
#include "XYRectangle.hpp"
#include "InstructionQueue.hpp"

#include "FeasibleNewton.hpp"
#include "TrustRegion.hpp"
#include "QuasiNewton.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
typedef Mesquite::SteepestDescent SolverType;

#include "TerminationCriterion.hpp"
#include "PMeanPTemplate.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "ElementPMeanP.hpp"
#include "MsqTimer.hpp"
#include "LaplacianSmoother.hpp"

#include <iostream>

using namespace Mesquite;

const double EPS = 1e-3;

bool find_sphere( Mesh* mesh, SphericalDomain* dom, MsqError& err );
bool find_cyl( Mesh* mesh, CylinderDomain* dom, MsqError& err );
bool find_cone( Mesh* mesh, ConicDomain* dom, MsqError& err );
void validate_hessian( Mesh* mesh, MeshDomain *geom );
 
class FDMetric : public TMPQualityMetric 
{
  private:
    bool fdGrad;
  public:
    FDMetric(TargetCalculator* tc, TargetMetric2D* mu, bool numerical_grad = true)
      : TMPQualityMetric( tc, mu, 0 ), fdGrad(numerical_grad) {}
    
    bool evaluate_with_gradient( PatchData& pd,
                                 size_t handle,
                                 double& value,
                                 std::vector<size_t>& indices,
                                 std::vector<Vector3D>& gradient,
                                 MsqError& err )
    {
      if (fdGrad)
        return QualityMetric::evaluate_with_gradient( pd, handle, value, indices, gradient, err );
      else
        return TMPQualityMetric::evaluate_with_gradient( pd, handle, value, indices, gradient, err );
    }
    
    bool evaluate_with_Hessian_diagonal( PatchData& pd,
                                         size_t handle,
                                         double& value,
                                         std::vector<size_t>& indices,
                                         std::vector<Vector3D>& gradient,
                                         std::vector<SymMatrix3D>& Hessian_diagonal,
                                         MsqError& err )
    {
      return QualityMetric::evaluate_with_Hessian_diagonal( pd, handle,
                               value, indices, gradient, Hessian_diagonal, err );
    }
    
    bool evaluate_with_Hessian( PatchData& pd,
                                size_t handle,
                                double& value,
                                std::vector<size_t>& indices,
                                std::vector<Vector3D>& gradient,
                                std::vector<Matrix3D>& Hessian,
                                MsqError& err )
    { return QualityMetric::evaluate_with_Hessian( pd, handle, value, indices,
                                                   gradient, Hessian, err );
    }
};
  
class FDTargetMetric : public TMetric 
{
  private:
    bool fdGrad;
  public:
    FDTargetMetric(bool numerical_grad = true) : fdGrad(numerical_grad) {}
    
    bool evaluate_with_grad( const MsqMatrix<2,2>& A,
                             const MsqMatrix<2,2>& W,
                             double& result,
                             MsqMatrix<2,2>& deriv_wrt_A,
                             MsqError& err )
    {
      if (fdGrad)
        return TargetMetric2D::evaluate_with_grad( A, W, result, deriv_wrt_A, err );
      else
        return TMetric::evaluate_with_grad( A, W, result, deriv_wrt_A, err );
    }
    
    bool evaluate_with_hess( const MsqMatrix<2,2>& A,
                    const MsqMatrix<2,2>& W,
                    double& result,
                    MsqMatrix<2,2>& deriv_wrt_A,
                    MsqMatrix<2,2> second_wrt_A[3],
                    MsqError& err )
    { return TargetMetric2D::evaluate_with_hess( A, W, result, deriv_wrt_A, second_wrt_A, err );
    }
};


void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " <input_file> [-o <output_file>] [-i <init_file>] [-p <plot_file>] [-t <timestep_base>] [-L <laplace_output>] " << std::endl;
  exit(1);
}

int main( int argc, char* argv[] )
{
  const char *input_file = 0, *output_file = 0, *init_file = 0, 
             *plot_file = 0, *timestep_base = 0, *laplace_file = 0;
  const char** arg = 0;
  for (int i = 1; i < argc; ++i) {
    if (arg) {
      if (*arg) {
        std::cerr << "File specified multiple times" <<std::endl;
        usage(argv[0]);
      }
      *arg = argv[i];
      arg = 0;
    } 
    else if (argv[i][0] == '-') {
      if (!argv[i][1] || argv[i][2])
        usage(argv[0]);
      switch (argv[i][1]) {
        case 'o': arg = &output_file;   break;
        case 'i': arg = &init_file;     break;
        case 'p': arg = &plot_file;     break;
        case 't': arg = &timestep_base; break;
        case 'L': arg = &laplace_file;  break;
        default: usage(argv[0]);
      }
    }
    else if (!input_file)
      input_file = argv[i];
    else 
      usage(argv[0]);
  }
  if (arg || !input_file) usage(argv[0]);

  MsqPrintError err(std::cerr);
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (err) {
    std::cerr << input_file << " : failed to read input file" << std::endl;
    return 2;
  }

  SphericalDomain sphere( Vector3D(0,0,0), 0.0 );
  CylinderDomain cylinder( 0, Vector3D(1,0,0), Vector3D(0,0,0) );
  ConicDomain cone;
  MeshDomain* geom = 0;
  if (find_sphere(&mesh, &sphere, err)) {
    geom = &sphere;
    std::cout << "Using spherical domain" << std::endl;
  }
  else if (err) 
    return 3;
  else if (find_cyl(&mesh, &cylinder, err)) {
    geom = &cylinder;
    std::cout << "Using cylinderical domain" << std::endl;
    mesh.mark_skin_fixed( err );
    if (err) {
      std::cerr << "Skinning failed" << std::endl;
      return 3;
    }
  }
  else if (err)
    return 3;
  else if (find_cone(&mesh, &cone, err)) {
    geom = &cone;
    std::cout << "Using conic domain" << std::endl;
    mesh.mark_skin_fixed( err );
    if (err) {
      std::cerr << "Skinning failed" << std::endl;
      return 3;
    }
  }
  else if (err)
    return 3;
  else {
    std::cerr << "Failed to identify mesh domain" << std::endl;
    return 2;
  }

  IdealShapeTarget tc;

//#ifdef NUMERICAL_HESSIAN
//  FDTargetMetric mu(false);
//#elif defined(NUMERICAL_GRADIENT)
//  FDTargetMetric mu(true);
//#else
  TMetric mu;
//#endif

#ifdef NUMERICAL_HESSIAN
  FDMetric qm( &tc, &mu, false );
#elif defined(NUMERICAL_GRADIENT)
  FDMetric qm( &tc, &mu, true );
#else
  TMPQualityMetric qm( &tc, &mu, 0 );
#endif
  ElementPMeanP avg_qm( 1.0, &qm );

  //validate_hessian( &mesh, geom );

  PMeanPTemplate of(1.0, &qm);
  SolverType qi(&of);
  TerminationCriterion inner;
  inner.add_absolute_vertex_movement( 1e-4 );
  if (plot_file) 
    inner.write_iterations( plot_file, err );
  if (timestep_base)
    inner.write_mesh_steps( timestep_base );
  
  qi.set_inner_termination_criterion(&inner);
  IdealWeightInverseMeanRatio imr;
  QualityAssessor qa(&qm, 0, 0.0, false, TMetric_NAME "_init", true, "Inverted_init");
  QualityAssessor qa2(&qm, 0, 0.0, false, TMetric_NAME, true, "Inverted");
  qa.add_quality_assessment(&avg_qm, 0, 0, "Average " TMetric_NAME "_init");
  qa2.add_quality_assessment(&avg_qm, 0, 0, "Average " TMetric_NAME );
  qa.add_quality_assessment(&imr, 0, 0, "InverseMeanRatio_init");
  qa2.add_quality_assessment(&imr, 0, 0, "InverseMeanRatio");
  InstructionQueue q1;
  q1.add_quality_assessor(&qa,err);
  //q1.run_instructions( &mesh, geom, err );
  if (init_file) {
    q1.run_instructions( &mesh, geom, err );
    if (err) {
      std::cerr << "Assessment failed!!!" << std::endl;
      return 3;
    }
    mesh.write_vtk( init_file, err );
    if (err) {
      std::cerr << init_file << " : failed to write output file" << std::endl;
      return 2;
    }
    std::cout << "Wrote initial mesh to file: " << init_file << std::endl;
  }

  InstructionQueue q2;
  q2.set_master_quality_improver(&qi,err);
  Timer timer;
  q2.run_instructions( &mesh, geom, err );
  if (err) {
    std::cerr << "Optimization failed!!!" << std::endl;
    return 3;
  }
  std::cout << "InstructionQueue returned after " << timer.since_birth() << " seconds" << std::endl;
  
  InstructionQueue q3;
  q3.add_quality_assessor(&qa2,err);
  q3.run_instructions( &mesh, geom, err );
  
  if (output_file) {
    mesh.write_vtk( output_file, err );
    if (err) {
      std::cerr << output_file << " : failed to write output file" << std::endl;
      return 2;
    }
    std::cout << "Wrote final mesh to file: " << output_file << std::endl;
  }
  
  if (laplace_file) {
    MeshImpl mesh2;
    mesh2.read_vtk( input_file, err );
    if (err) {
      std::cerr << input_file << " : failed to read input file" << std::endl;
      return 2;
    }
    
    if (geom != &sphere) {
      mesh2.mark_skin_fixed( err );
      if (err) {
        std::cerr << "Skinning failed" << std::endl;
        return 3;
      }
    }
    
    std::cout << std::endl << "Laplace smooth:" << std::endl;
    
    InstructionQueue q4;
    LaplacianSmoother laplace;
    laplace.set_outer_termination_criterion(&inner);
    q4.set_master_quality_improver( &laplace, err );
    q4.add_quality_assessor(&qa2,err);
    q4.run_instructions( &mesh2, geom, err );
    if (err) {
      std::cerr << "Laplace smooth failed!!!" << std::endl;
      return 3;
    }
    
    mesh2.write_vtk( laplace_file, err );
    if (err) {
      std::cerr << laplace_file << " : failed to write output file" << std::endl;
      return 2;
    }
    std::cout << "Wrote Laplace smoothed mesh to file: " << laplace_file << std::endl;
  }
  
  return 0;
}

bool find_sphere( Mesh* mesh, SphericalDomain* dom, MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  Vector3D min(coords[0]), max(coords[0]);
  for (size_t i = 1; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
  
  max *= 0.5;
  min *= 0.5;
  for (int j = 0; j < 3; ++j) {
    max[j] = round(20*max[j])/20;
    min[j] = round(20*min[j])/20;
  }
  Vector3D diff = max - min;
  Vector3D cent = max + min;
  double rad = diff[0];
  bool ok = true;
  for (size_t i = 0; i < coords.size(); ++i) {
    diff = coords[i]-cent;
    if (fabs(diff.length() - rad) > EPS) {
      //std::cerr << "Vertex " << i << " is not on mean sphere" << std::endl;
      //std::cerr << "  cent: " << cent << ", rad: " << rad << ", coords: " << coords[i] << std::endl;
      ok = false;
    }
  }

 dom->set_sphere( cent, rad );
 return ok;
}
  
bool find_cyl( Mesh* mesh, CylinderDomain* dom, MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  Vector3D min(coords[0]), max(coords[0]);
  for (size_t i = 1; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
  
  Vector3D mid = 0.5*(min+max);
  min.set(HUGE_VAL,HUGE_VAL,HUGE_VAL);
  max.set(0,0,0);
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D diff = mid - coords[i];
    for (int j = 0; j < 3; ++j) {
      double x1 = diff[(j+1)%3];
      double x2 = diff[(j+2)%3];
      double dist_sqr = x1*x1 + x2*x2;
      if (min[j] > dist_sqr)
        min[j] = dist_sqr;
      if (max[j] < dist_sqr)
        max[j] = dist_sqr;
    }
  }
  
  max -= min;
  int axis = (max[0] < max[1]) && (max[0] < max[2]) ? 0 : (max[1] < max[2]) ? 1 : 2;
  double rad = round(20*(sqrt(min[axis]+0.5*max[axis])))/20;
  for (int j = 0; j < 3; ++j) 
    mid[j] = round(20*mid[j])/20;
  
  bool ok = true;
  const int a = (axis+1)%3, b = (axis+2)%3;
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D diff = mid - coords[i];
    for (int j = 0; j < 3; ++j) {
      double dist = sqrt(diff[a]*diff[a]+diff[b]*diff[b]);
      if (fabs(dist - rad) > EPS)
        ok = false;
    }
  }
  
  Vector3D axis_vec(0,0,0);
  axis_vec[axis] = 1.0;
  *dom = CylinderDomain( rad, axis_vec, mid );
  return ok;
}
  
bool find_cone( Mesh* mesh, ConicDomain* dom, MsqError& err )
{
  ConicDomain cone( 2, 10 );

  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D p(coords[i]);
    cone.snap_to( 0, p );
    if ((p - coords[i]).length_squared() > 1e-5)
      return false;
  }
  
  *dom = cone;
  return true;
}

void validate_hessian( Mesh* mesh, MeshDomain *geom )
{
  

  IdealShapeTarget tc;
  TMetric mu;
  std::vector<size_t> indices1, indices2;
  std::vector<Vector3D> grad1, grad2;
  std::vector<Matrix3D> hess1, hess2;
  double val1, val2;

#ifdef NUMERICAL_HESSIAN
  FDMetric nqm( &tc, &mu, false );
#else
  FDMetric nqm( &tc, &mu, true );
#endif
  TMPQualityMetric qm( &tc, &mu, 0 );

  MsqPrintError err(std::cerr);
  PatchData pd;
  pd.set_mesh( mesh );
  pd.set_domain( geom );
  pd.fill_global_patch( err );
  
  std::vector<size_t>::iterator j;
  std::vector<size_t> metric_handles;
  qm.get_evaluations( pd, metric_handles, false, err );
  for (j = metric_handles.begin(); j != metric_handles.end(); ++j) 
  {
    indices1.clear(); grad1.clear(); hess1.clear();
    qm.evaluate_with_Hessian( pd, *j, val1, indices1, grad1, hess1, err );
    assert(!err);
    indices2.clear(); grad2.clear(); hess2.clear();
    nqm.evaluate_with_Hessian( pd, *j, val2, indices2, grad2, hess2, err );
    assert(!err);
    assert(indices1 == indices2);
    assert(fabs(val1 - val2) < EPS);
    bool OK = true;
    for (size_t i = 0; i < grad1.size(); ++i)
      for (size_t d = 0; d < 3; ++d)
        if (fabs(grad1[i][d] - grad2[i][d]) > 100*EPS) {
          std::cerr << "Gradient term for index " << indices1[i] << " axis " << d
                    << " for sample " << *j << " has analytic value of " 
                    << grad1[i][d] << " and numerical value of "
                    << grad2[i][d] << std::endl;
          OK = false;
        }
    size_t h = 0;
    for (size_t i = 0; i < indices1.size(); ++i)
      for (size_t k = i; k < indices1.size(); ++k, ++h)
        for (size_t r = 0; r < 3; ++r)
          for (size_t c = 0; c < 3; ++c)
            if (fabs(hess1[h][r][c] - hess2[h][r][c]) > 100*EPS) {
              std::cerr << "Hessian term " << r << "," << c << " for indices " 
                        << indices1[i] << ", " << indices1[k] 
                        << " for sample " << *j 
                        << " has analytic value of " << hess1[h][r][c] 
                        << " and numerical value of " << hess2[h][r][c] << std::endl;
              OK = false;
            }
      
    if (!OK) {
      qm.evaluate_with_Hessian( pd, *j, val1, indices1, grad1, hess1, err );
    } 
  }
}
