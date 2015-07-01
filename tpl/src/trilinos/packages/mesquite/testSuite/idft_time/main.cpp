/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    kraftche@cae.wisc.edu     
   
  ***************************************************************** */

#include "meshfiles.h"

#include <iostream>
#include <sstream>
using std::cout;
using std::cerr;
using std::endl;

#include <stdlib.h>

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"

#include "IdealWeightInverseMeanRatio.hpp"
#include "ElementPMeanP.hpp"
#include "TInverseMeanRatio.hpp"
#include "TQualityMetric.hpp"
#include "TOffset.hpp"
#include "TOffset.hpp"

#include "IdealShapeTarget.hpp"
#include "PMeanPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "QuasiNewton.hpp"

#include "MeshWriter.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;

const char* DEFAULT_INPUT = MESH_FILES_DIR "2D/vtk/tris/untangled/equil_tri2.vtk";

/* Print usage or help information: exits if err == true */
void usage( const char* argv0 = 0, bool err = true )
{
  std::ostream& s = err ? cerr : cout;
  const char* defname = "main";
  if (!argv0)
    argv0 = defname;
  
  s << "Usage: " << defname << " [-n|-c|-q] [-e] [-t] [-a] [-N] [{-v|-g|-p} <output_file>] [<input_file>]" << endl
    << "       " << defname << " -h" << std::endl;
  
  if (err)
    exit(1);

  s << "  -n : Use SteepestDescent solver (default)" << endl
    << "  -c : Use ConjugateGradient solver" << endl
    << "  -q : Use QuasiNewton solver" << endl
    << "  -e : Test IdealWeightInverseMeanRatio metric" << endl
    << "  -t : Test InverseMeanRatio3D target metric" << endl
    << "  -a : Test ElementPMeanP(InverseMeanRatio3D)" << endl
    << "  -N : Test InverseMeanRatio3D with finite difference derivatives" << endl
    << "  -C : Compare InverseMeanRatio3D and IdealWeightInverseMeanRatio" << endl
    << "  -v : Write output mesh to specified VTK file" << endl
    << "  -g : Write output mesh to specified GnuPlot file" << endl
    << "  -p : Write output mesh to specified EPS file" << endl
    << "  -P : Write solver plot data to specified file." << endl
    << "Default is all metrics if non specified." << endl
    << "Default input file: " << DEFAULT_INPUT << endl;
}

const char* vtk_file = 0; /* vtk output file name */
const char* eps_file = 0; /* eps output file name */
const char* gpt_file = 0; /* GNUPlot output file name */
const char* plot_file = 0; /* Time-dependent plot of solver data */

enum Solver { STEEP_DESCENT, CONJ_GRAD, QUASI_NEWT };

/* Run an optimization: returns average quality of final mesh */
double run( QualityMetric* metric, 
            Solver solver,
            const char* input_file,
            double& seconds_out,
            int& iterations_out );

/* Force finite difference approximation of derivatives for 
 * an existing quality metric */
class NumericQM : public QualityMetric {
  private:
    QualityMetric* realMetric;
  public:
    NumericQM( QualityMetric* real_metric ) : realMetric(real_metric) {}
    virtual MetricType get_metric_type() const;
    virtual std::string get_name() const;
    virtual int get_negate_flag() const;
    virtual void get_evaluations( PatchData& pd,
                          std::vector<size_t>& handles, 
                          bool free_vertices_only,
                          MsqError& err );
    virtual bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );
    virtual bool evaluate_with_indices( PatchData& pd,
                   size_t handle,
                   double& value,
                   std::vector<size_t>& indices,
                   MsqError& err );
};
QualityMetric::MetricType NumericQM::get_metric_type() const
  { return realMetric->get_metric_type(); }
std::string NumericQM::get_name() const{
  std::string r = realMetric->get_name();
  r += " (FD)";
  return r;
}

/* At each evaluation of this metric, compare the values resulting
 * from the evaluation of two other metrics: flag an error if they
 * differ or return the common result if the same */
class CompareMetric : public QualityMetric 
{
private:
  QualityMetric *metric1, *metric2;
  bool maskPlane;
  int maskAxis;
  std::vector<size_t> m2Handles;
  std::vector<Vector3D> m2Grad;
  std::vector<SymMatrix3D> m2Diag;
  std::vector<Matrix3D> m2Hess;
  static const double epsilon;
public:
  CompareMetric( QualityMetric* qm1, 
                 QualityMetric* qm2,
                 bool mask_qm2_coord = false )
    : metric1(qm1), metric2(qm2), maskPlane(mask_qm2_coord), maskAxis(-1)
    {}
  
  MetricType get_metric_type() const;
  
  std::string get_name() const;
  
  int get_negate_flag() const;
  
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
                        
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );
                 
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err );
                 
  bool evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err );
                 
  bool evaluate_with_Hessian_diagonal( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<SymMatrix3D>& Hessian_diagonal,
                 MsqError& err );
                 
  bool evaluate_with_Hessian( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<Matrix3D>& Hessian,
                 MsqError& err );
                 
  void get_mask_axis( PatchData& pd );
  bool equal(    Vector3D grad1, const    Vector3D& grad2 ) const;
  bool equal( SymMatrix3D hess1, const SymMatrix3D& hess2 ) const;
  bool equal(    Matrix3D hess1, const    Matrix3D& hess2 ) const;
};
const double CompareMetric::epsilon = 5e-2;
  
/* Parse command line options and call 'run' */
int main( int argc, char* argv[] )
{
  Solver solver = STEEP_DESCENT;
  bool do_non_target_metric = false;
  bool do_new_target_metric = false;
  bool do_new_target_average = false;
  bool do_new_target_numeric = false;
  bool do_compare_metric = false;
  const char* input_file = DEFAULT_INPUT;
  bool no_more_flags = false;
  
  std::list<const char**> exp_list;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-' && !no_more_flags) {
      for (int k = 1; argv[i][k]; ++k) {
        switch (argv[i][k]) {
          case 'n': solver = STEEP_DESCENT; break;
          case 'c': solver = CONJ_GRAD; break;
          case 'q': solver = QUASI_NEWT; break;
          case 'e': do_non_target_metric = true; break;
          case 't': do_new_target_metric = true; break;
          case 'a': do_new_target_average = true; break;
          case 'N': do_new_target_numeric = true; break;
          case 'C': do_compare_metric = true; break;
          case 'v': exp_list.push_back(&vtk_file); break;
          case 'g': exp_list.push_back(&gpt_file); break;
          case 'p': exp_list.push_back(&eps_file); break;
          case 'P': exp_list.push_back(&plot_file); break;
          case 'h': usage(argv[0],false); return 0;
          case '-': no_more_flags = true; break;
          default:
            cerr << "Invalid flag: '" << argv[i][k] << "'" << endl;
            usage(argv[0]);
        }
      }
    }
    else if (!exp_list.empty()) {
      const char** ptr = exp_list.front();
      exp_list.pop_front();
      *ptr = argv[i];
    }
    else {
      if (input_file != DEFAULT_INPUT) {
        cerr << "Unexpected argument: \"" << argv[i] << '"' << endl;
        usage(argv[0]);
      }
      input_file = argv[i];
    }
  }
  
  int count = 0;
  if (do_non_target_metric)
    ++count;
  if (do_new_target_metric)
    ++count;
  if (do_new_target_average)
    ++count;
  if (do_new_target_numeric)
    ++count;
  if (do_compare_metric)
    ++count;
  
  if (!count) {
    do_compare_metric = true;
    count = 1;
  }
  
  if ((vtk_file || gpt_file || eps_file) && count != 1) {
    cerr << "Error: Cannot write output file if running multiple tests" << endl;
    return 2;
  }
  
  IdealWeightInverseMeanRatio non_target_metric;
  IdealShapeTarget new_target;
  TInverseMeanRatio tmp;
  TOffset tmp_off( 1.0, &tmp ); // target metrics are zero-based
  TQualityMetric new_target_metric( &new_target, &tmp_off );
  ElementPMeanP new_target_average( 1.0, &new_target_metric );
  NumericQM new_target_numeric( &new_target_metric );
  CompareMetric comp_metric( &non_target_metric, &new_target_average, true );
  
  std::ostringstream os;
  double secs,qual;
  if (do_non_target_metric) {
    qual = run( &non_target_metric, solver, input_file, secs, count );
    os << "IdealWeightInverseMeanRatio: " << qual << " after " << count << " iterations in " << secs << " seconds" << endl;
  }
  if (do_new_target_metric) {
    qual = run( &new_target_metric, solver, input_file, secs, count );
    os << "TQualityMetric              : " << qual << " after " << count << " iterations in " << secs << " seconds" << endl;
  }
  if (do_new_target_average) {
    qual = run( &new_target_average, solver, input_file, secs, count );
    os << "ElementPMeanP              : " << qual << " after " << count << " iterations in " << secs << " seconds" << endl;
  }
  if (do_new_target_numeric) {
    qual = run( &new_target_numeric, solver, input_file, secs, count );
    os << "TQualityMetric (FD)         : " << qual << " after " << count << " iterations in " << secs << " seconds" << endl;
  }
  if (do_compare_metric) {
    qual = run( &comp_metric, solver, input_file, secs, count );
    os << "Metric comparison      : " << qual << " after " << count << " iterations in " << secs << " seconds" << endl;
  }
  
  cout << endl << os.str() << endl;
  return 0;
}



double run( QualityMetric* metric, 
            Solver solver_type,
            const char* input_file,
            double& seconds_out,
            int& iterations_out )
{
  MsqPrintError err(cerr);
  IdealWeightInverseMeanRatio qa_metric;
  TerminationCriterion inner, outer;
  outer.add_iteration_limit( 1 );
  inner.add_absolute_vertex_movement( 1e-4 );
  inner.add_iteration_limit( 100 );
  PMeanPTemplate of( 1.0, metric );
  QualityAssessor qa( &qa_metric );
  qa.add_quality_assessment( metric );
  InstructionQueue q;
  SteepestDescent steep(&of);
  QuasiNewton quasi(&of);
  ConjugateGradient conj(&of);
  VertexMover* solver = 0;
  switch (solver_type) {
    case STEEP_DESCENT: solver = &steep; break;
    case QUASI_NEWT:solver = &quasi;break;
    case CONJ_GRAD: solver = &conj; break;
  }
  q.set_master_quality_improver( solver, err );
  q.add_quality_assessor( &qa, err );
  solver->set_inner_termination_criterion(&inner);
  solver->set_outer_termination_criterion(&outer);
  
  
  if (plot_file)
    inner.write_iterations( plot_file, err );
  
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (err) {
    cerr << "Failed to read input file: \"" << input_file << '"' << endl;
    exit(1);
  }
  
  std::vector<Mesh::VertexHandle> handles;
  mesh.get_all_vertices( handles, err );
  if (handles.empty()) {
    cerr << "no veritces in mesh" << endl;
    exit(1);
  }
  std::vector<MsqVertex> coords(handles.size());
  mesh.vertices_get_coordinates( arrptr(handles), arrptr(coords), handles.size(), err );
  Vector3D min(HUGE_VAL), max(-HUGE_VAL);
  for (size_t i = 0; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
  
  Vector3D size = max - min;
  PlanarDomain* domain = 0;
  if (size[0] < 1e-4) 
    domain = new PlanarDomain( PlanarDomain::YZ, min[0] );
  else if (size[1] < 1e-4)
    domain = new PlanarDomain( PlanarDomain::XZ, min[1] );
  else if (size[2] < 1e-4)
    domain = new PlanarDomain( PlanarDomain::XY, min[2] );
  
  Timer timer;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
  q.run_instructions( &mesh_and_domain, err );
  seconds_out = timer.since_birth();
  if (err) {
    cerr << "Optimization failed." << endl << err << endl;
    abort();
  }

  if (vtk_file) {
    MeshWriter::write_vtk( &mesh, vtk_file, err );
    if (err) 
      cerr << vtk_file << ": failed to write file." << endl;
  }
  if (gpt_file) {
    MeshWriter::write_gnuplot( &mesh, gpt_file, err );
    if (err) 
      cerr << gpt_file << ": failed to write file." << endl;
  }
  if (eps_file) {
    PlanarDomain xy(PlanarDomain::XY);
    MeshWriter::Projection proj( domain ? domain : &xy );
    MeshWriter::write_eps( &mesh, eps_file, proj, err );
    if (err) 
      cerr << eps_file << ": failed to write file." << endl;
  }
  delete domain;
  
  iterations_out = inner.get_iteration_count();
  
  const QualityAssessor::Assessor* a = qa.get_results( &qa_metric );
  return a->get_average();
}
  
int NumericQM::get_negate_flag() const
{ return realMetric->get_negate_flag(); }

void NumericQM::get_evaluations( PatchData& pd,
                      std::vector<size_t>& handles, 
                      bool free_vertices_only,
                      MsqError& err )
{ realMetric->get_evaluations( pd, handles, free_vertices_only, err ); }

bool NumericQM::evaluate( PatchData& pd, 
                size_t handle, 
                double& value, 
                MsqError& err )
{ return realMetric->evaluate( pd, handle, value, err ); }

bool NumericQM::evaluate_with_indices( PatchData& pd,
               size_t handle,
               double& value,
               std::vector<size_t>& indices,
               MsqError& err )
{ return realMetric->evaluate_with_indices( pd, handle, value, indices, err ); }


QualityMetric::MetricType CompareMetric::get_metric_type() const
{
  MetricType t1 = metric1->get_metric_type();
  assert(metric2->get_metric_type() == t1);
  return t1;
}
  
std::string CompareMetric::get_name() const
{
  std::string n = metric1->get_name();
  n += " =? ";
  n += metric2->get_name();
  return n;
}
  
int CompareMetric::get_negate_flag() const
{
  assert(metric1->get_negate_flag() == metric2->get_negate_flag());
  return metric1->get_negate_flag();
}
  
void CompareMetric::get_evaluations( PatchData& pd, 
                                     std::vector<size_t>& handles, 
                                     bool free_vertices_only,
                                     MsqError& err )
{
  if (maskPlane)
    get_mask_axis(pd);

  m2Handles.clear();
  metric1->get_evaluations( pd, handles, free_vertices_only, err ); MSQ_ERRRTN(err);
  metric2->get_evaluations( pd, m2Handles, free_vertices_only, err ); MSQ_ERRRTN(err);
  bool same = (handles.size() == m2Handles.size());
  std::sort( m2Handles.begin(), m2Handles.end() );
  for (std::vector<size_t>::iterator i = handles.begin(); i != handles.end(); ++i) 
    if (!std::binary_search( m2Handles.begin(), m2Handles.end(), *i ))
      same = false;
  if (!same) {
    MSQ_SETERR(err)("Metrics have incompatible lists of evaluation handles.\n",
                    MsqError::INVALID_STATE);
  }
}


bool CompareMetric::evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err )
{
  double m2val;
  bool r1, r2;
  r1 = metric1->evaluate( pd, handle, value, err ); MSQ_ERRZERO(err);
  r2 = metric2->evaluate( pd, handle, m2val, err ); MSQ_ERRZERO(err);
  if (r1 != r2 || (r1 && fabs(value - m2val) > epsilon)) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
                    "Metrics returned different values for handle %lu "
                    "in evaluate:\n"
                    "\t%s %f vs. %s %f\n", (unsigned long)handle,
                    r1?"true":"false",value,r2?"true":"false",m2val);
  }
  
  return r1 && !err;
}


bool CompareMetric::evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err )
{
  double m2val;
  bool r1, r2;
  m2Handles.clear();
  r1 = metric1->evaluate_with_indices( pd, handle, value, indices, err ); MSQ_ERRZERO(err);
  r2 = metric2->evaluate_with_indices( pd, handle, m2val, m2Handles, err ); MSQ_ERRZERO(err);
  if (r1 != r2 || (r1 && fabs(value - m2val) > epsilon)) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
                    "Metrics returned different values for handle %lu "
                    "in evaluate_with_indices:\n"
                    "\t%s %f vs. %s %f\n", (unsigned long)handle,
                    r1?"true":"false",value,r2?"true":"false",m2val);
  }
  else {
    bool same = (indices.size() == m2Handles.size());
    std::sort( m2Handles.begin(), m2Handles.end() );
    for (std::vector<size_t>::iterator i = indices.begin(); i != indices.end(); ++i) 
      if (!std::binary_search( m2Handles.begin(), m2Handles.end(), *i ))
        same = false;
    if (!same) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned incompatible lists of vertex indices"
                      " for handle %lu in evaluate_with_indices\n.", 
                      (unsigned long)handle );
    }
  }
  
  return r1 && !err;
}


bool CompareMetric::evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err )
{
  double m2val;
  bool r1, r2;
  m2Handles.clear();
  m2Grad.clear();
  r1 = metric1->evaluate_with_gradient( pd, handle, value, indices, gradient, err ); MSQ_ERRZERO(err);
  r2 = metric2->evaluate_with_gradient( pd, handle, m2val, m2Handles, m2Grad, err ); MSQ_ERRZERO(err);
  if (r1 != r2 || (r1 && fabs(value - m2val) > epsilon)) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
                    "Metrics returned different values for handle %lu in "
                    "evaluate_with_gradient:\n"
                    "\t%s %f vs. %s %f\n", (unsigned long)handle,
                    r1?"true":"false",value,r2?"true":"false",m2val);
  }
  else {
    std::vector<size_t>::const_iterator i, j;
    std::vector<Vector3D>::const_iterator r, s;
    int grad_diff = 0;
    bool same = (indices.size() == m2Handles.size());
    std::sort( m2Handles.begin(), m2Handles.end() );
    for (i = indices.begin(); i != indices.end(); ++i) {
      j = std::lower_bound( m2Handles.begin(), m2Handles.end(), *i );
      if (j == m2Handles.end() || *j != *i) {
        same = false;
        continue;
      }
      
      r = gradient.begin() + (i - indices.begin());
      s = m2Grad.begin() + (j - m2Handles.begin());
      if (!equal(*r,*s))
        ++grad_diff;
    }
      
    if (!same) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned incompatible lists of vertex indices"
                      " for handle %lu in evaluate_with_gradient\n.", 
                      (unsigned long)handle );
    }
    else if (grad_diff) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned different gradient vectors for "
                      " %d of %u vertices for handle %lu in "
                      "evaluate_with_gradient\n.", 
                      grad_diff, (unsigned)gradient.size(), 
                      (unsigned long)handle );
    }
  }
  
  return r1 && !err;
}

bool CompareMetric::evaluate_with_Hessian_diagonal( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<SymMatrix3D>& diagonal,
                 MsqError& err )
{
  double m2val;
  bool r1, r2;
  m2Handles.clear();
  m2Grad.clear();
  m2Diag.clear();
  r1 = metric1->evaluate_with_Hessian_diagonal( pd, handle, value, indices, gradient, diagonal, err ); MSQ_ERRZERO(err);
  r2 = metric2->evaluate_with_Hessian_diagonal( pd, handle, m2val, m2Handles, m2Grad, m2Diag, err ); MSQ_ERRZERO(err);
  if (r1 != r2 || (r1 && fabs(value - m2val) > epsilon)) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
                    "Metrics returned different values for handle %lu in "
                    "evaluate_with_Hessian_diagonal:\n"
                    "\t%s %f vs. %s %f\n", (unsigned long)handle,
                    r1?"true":"false",value,r2?"true":"false",m2val);
  }
  else {
    std::vector<size_t>::const_iterator i, j;
    std::vector<Vector3D>::const_iterator r, s;
    std::vector<SymMatrix3D>::const_iterator u, v;
    int grad_diff = 0, hess_diff = 0;
    bool same = (indices.size() == m2Handles.size());
    std::sort( m2Handles.begin(), m2Handles.end() );
    for (i = indices.begin(); i != indices.end(); ++i) {
      j = std::lower_bound( m2Handles.begin(), m2Handles.end(), *i );
      if (j == m2Handles.end() || *j != *i) {
        same = false;
        continue;
      }
      
      r = gradient.begin() + (i - indices.begin());
      s = m2Grad.begin() + (j - m2Handles.begin());
      if (!equal(*r,*s))
        ++grad_diff;
      
      u = diagonal.begin() + (i - indices.begin());
      v = m2Diag.begin() + (j - m2Handles.begin());
      if (!equal(*u,*v))
        ++hess_diff;
    }
      
    if (!same) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned incompatible lists of vertex indices"
                      " for handle %lu in evaluate_with_Hessian_diagonal\n.", 
                      (unsigned long)handle );
    }
    else if (grad_diff) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned different gradient vectors for "
                      " %d of %u vertices for handle %lu in "
                      "evaluate_with_Hessian_diagonal\n.", 
                      grad_diff, (unsigned)gradient.size(), 
                      (unsigned long)handle );
    }
    else if (hess_diff) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned different Hessian blocks for "
                      " %d of %u vertices for handle %lu in "
                      "evaluate_with_Hessian_diagonal\n.", 
                      hess_diff, (unsigned)diagonal.size(), 
                      (unsigned long)handle );
    }
  }
  
  return r1 && !err;
}
                 
bool CompareMetric::evaluate_with_Hessian( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<Matrix3D>& Hessian,
                 MsqError& err )
{
  double m2val;
  bool r1, r2;
  m2Handles.clear();
  m2Grad.clear();
  m2Hess.clear();
  r1 = metric1->evaluate_with_Hessian( pd, handle, value, indices, gradient, Hessian, err ); MSQ_ERRZERO(err);
  r2 = metric2->evaluate_with_Hessian( pd, handle, m2val, m2Handles, m2Grad, m2Hess, err ); MSQ_ERRZERO(err);
  if (r1 != r2 || fabs(value - m2val) > epsilon) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
                    "Metrics returned different values for handle %lu in "
                    "evaluate_with_Hessian:\n"
                    "\t%s %f vs. %s %f\n", (unsigned long)handle,
                    r1?"true":"false",value,r2?"true":"false",m2val);
  }
  else {
    std::vector<size_t>::const_iterator i, j;
    std::vector<Vector3D>::const_iterator r, s;
    int grad_diff = 0, hess_diff = 0;
    bool same = (indices.size() == m2Handles.size());
    std::sort( m2Handles.begin(), m2Handles.end() );
    for (i = indices.begin(); i != indices.end(); ++i) {
      j = std::lower_bound( m2Handles.begin(), m2Handles.end(), *i );
      if (j == m2Handles.end() || *j != *i) {
        same = false;
        continue;
      }
      
      r = gradient.begin() + (i - indices.begin());
      s = m2Grad.begin() + (j - m2Handles.begin());
      if (!equal(*r,*s)) {
        ++grad_diff;
          // call again for so debugger can step into it after failure is found
        std::vector<size_t> i2;
        std::vector<Vector3D> g2;
        std::vector<Matrix3D> h2;
        metric2->evaluate_with_Hessian(pd, handle, m2val, i2, g2, h2, err );
      }
    }  
     
    if (!same) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned incompatible lists of vertex indices"
                      " for handle %lu in evaluate_with_Hessian\n.", 
                      (unsigned long)handle );
    }
    else if (grad_diff) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,
                      "Metrics returned different gradient vectors for "
                      " %d of %u vertices for handle %lu in "
                      "evaluate_with_Hessian\n.", 
                      grad_diff, (unsigned)gradient.size(), 
                      (unsigned long)handle );
    }
    else {
      size_t row, col, row2, col2, idx, idx2;
      for (row = idx = 0; row < indices.size(); ++row) {
        row2 = std::lower_bound( m2Handles.begin(), m2Handles.end(), indices[row] ) - m2Handles.begin();
        for (col = row; col < indices.size(); ++col, ++idx) {
          col2 = std::lower_bound( m2Handles.begin(), m2Handles.end(), indices[col] ) - m2Handles.begin();
          if (row2 <= col2) {
            idx2 = indices.size()*row2 - row2*(row2+1)/2 + col2;
            if (!equal(Hessian[idx], m2Hess[idx2]))
              ++hess_diff;
          }
          else {
            idx2 = indices.size()*col2 - col2*(col2+1)/2 + row2;
            if (!equal(Hessian[idx], transpose(m2Hess[idx2])))
              ++hess_diff;
          }
        }
      }

      if (hess_diff) {
        MSQ_SETERR(err)(MsqError::INVALID_STATE,
                        "Metrics returned different Hessian blocks for "
                        " %d of %u vertices for handle %lu in "
                        "evaluate_with_Hessian\n.", 
                        hess_diff, (unsigned)Hessian.size(), 
                        (unsigned long)handle );
      }
    }
  }
  
  return r1 && !err;
}

void CompareMetric::get_mask_axis( PatchData& pd )
{
  maskAxis = -1;
  PlanarDomain* dom = reinterpret_cast<PlanarDomain*>(pd.get_domain());
  int bits = 0;
  if (dom) {
    Vector3D n = dom->get_normal();
    for (int i = 0; i < 3; ++i)
      if (fabs(n[i]) < epsilon)
        bits |= (1 << i);
    switch (bits) {
      case 3: maskAxis = 2; break;
      case 5: maskAxis = 1; break;
      case 6: maskAxis = 0; break;
    }
  }
}

bool CompareMetric::equal( Vector3D grad1, const Vector3D& grad2 ) const
{
  if (maskAxis >= 0)
    grad1[maskAxis] = grad2[maskAxis];
  return (grad1 - grad2).length_squared() <= epsilon*epsilon;
}

bool CompareMetric::equal( SymMatrix3D hess1, const SymMatrix3D& hess2 ) const
{
  if (maskAxis >= 0) {
    for (unsigned i = 0; i < 3; ++i) {
      hess1(maskAxis,i) = hess2(maskAxis,i);
      hess1(i,maskAxis) = hess2(i,maskAxis);
    }
  }
  return Frobenius_2(hess1-hess2) <= epsilon*epsilon;
}

bool CompareMetric::equal( Matrix3D hess1, const Matrix3D& hess2 ) const
{
  if (maskAxis >= 0) {
    for (unsigned i = 0; i < 3; ++i) {
      hess1[maskAxis][i] = hess2[maskAxis][i];
      hess1[i][maskAxis] = hess2[i][maskAxis];
    }
  }
  return Frobenius_2(hess1-hess2) <= epsilon*epsilon;
}
