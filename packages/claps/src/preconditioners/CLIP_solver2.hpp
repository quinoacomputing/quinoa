//@HEADER
// ************************************************************************
//
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef CLIP_SOLVER2_HPP
#define CLIP_SOLVER2_HPP
#include <mpi.h>
#include <math.h>
#include <algorithm>
//#include "CLIP_graph.hpp"
//#include "CLIP_sub.hpp"
#include "CLOP_constraint.hpp"
#include "CRD_utils.hpp"
//#include "../include/sort_prototypes.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "EpetraExtCD_MatrixMatrix.hpp"
#include "CRS_serial.hpp"
#include "preconditioner_crd.hpp"
#include "sparse_lu2.hpp"
#if defined(AMG_CRD)
  #include "amg_solve.hpp"
#endif


//class CLIP_solver2 : public preconditioner_crd
class CLIP_solver2 : public preconditioner_crd
{
 public: // functions
  CLIP_solver2(CRS_serial* A_,              // problem data
	       const Epetra_Map* SubMap_,   // subdomain map
	       const Epetra_Map* OwnMap_,   // unique owner map
	       const double* clip_params_,  // array of solver parameters
	       const double* amg_params_);  // array of amg parameters
  ~CLIP_solver2();
  double norm2(double a[], int n);
  double dotprod(double a[], double b[], int n);
  void sum_vectors(double a[], int n, double a_sum[]);
  int initialize_solve(double u[], double r[]);
  void apply_preconditioner(const double r[], double z[]);
  void A_times_x(double* x, double* Ax);
 private: // functions
  void determine_components();
  void determine_components(int A1[], int A2[], int N, 
			    int* & comp1, int* & comp2, int & ncomp);
  void determine_dof_sets();
  void determine_corner_dofs();
  void share_corner_info();
  void modify_dof_sets();
  void factor_sub_matrices();
  void calculate_coarse();
  void allocate_vectors();
  void static_correction1();
  void static_correction2();
  void coarse_correction();
  void substructure_correction();
  void weight_scale(Epetra_MultiVector* Mvec);
  void A_times_x_sub(double* x, double* Ax);
  void A_times_x_sub(double* x, double* Ax, int* rows, int n);
  void get_matrix_diag(CRS_serial* AA, double* & diag);
  int determine_shell(int node);
  bool find_sub(int dof, int sub);
  void check_two_nodes(int ii, int* adof, bool* adof_flag,
		       int* anode, bool* anode_flag);
  void check_three_nodes(int ii, int* adof, bool* adof_flag,
			 int* anode, bool* anode_flag);
  void get_adof(int ii, int* aa, int nn, int* adof, int & nadof, 
		bool* adof_flag, int tflag=0);
  bool contains(int ii, int dof2);
  void get_anode(int* adof, int nadof, int* anode, 
		 bool* anode_flag, int & nanode);
  void set_corner_flag(int node);
  int find_node1(int* anode, int nanode);
  void find_farthest1(int* anode, int nanode, int node1, int & node2);
  void find_farthest2(int* anode, int nanode, int & node1, int & node2);
  void find_node1_node2(int* anode, int nanode, int & node1, int & node2);
  int find_max_area(int* anode, int nanode, int node1, int node2);
  void gen_matrix(CRS_serial* A, int na, int adof[], 
		  int* & rowbeg, int* & colidx, double* &vals);
  void zero_pointers();
  void spmat_datfile(const Epetra_CrsMatrix & A, char fname[], int opt);
  void gen_amg_data(int n, int* dofa, double* & amg_mat, int* & amg_A1,
		    int* & amg_A2, double* & amg_x, double* & amg_y,
		    double* & amg_z, int* & amg_nodebeg,
		    int* & amg_local_dof, int & amg_nnode);

 private: // variables
  CRS_serial *A;
  const Epetra_Map *SubMap, *OwnMap;
  const double *clip_params, *amg_params;
  const Epetra_Comm & Comm;

  int cdof_option, maxiter, atype, ndim, local_solver, prt_debug, prt_summary;
  int max_orthog, chk_sub_singularity, krylov_method, scale_option, sub_solver;
  int num_rigid_mode, coarse_solver;
  double solver_tol;

  // extra to work for now
  const Epetra_IntVector *NodalDofs;
  const Epetra_CrsMatrix *ASub;
  Epetra_Map *StandardMap;
  // current stuff
  Epetra_Import *Importer, *ImporterC, *ImporterB;
  Epetra_Export *Exporter, *ExporterC, *ExporterB;
  Epetra_Map *SubMapC, *SubMapB;
  Epetra_Map *OwnMapC, *OwnMapB;
  Epetra_MultiVector *rOwn, *rSub, *rOwnB, *rSubB, *rOwnC, *rSubC;
  Epetra_MultiVector *uOwn, *uSub, *uOwnB, *uSubB, *uOwnC, *uSubC;
  int MyPID, NumProc, ndof_sub, ndof_global, print_flag;
  int *comp1, *comp2, ncomp, ncomp_sum, *sub1, *sub2, *dset1, *dset2;
  int ndof_set, *corner_flag, *dof2node, *local_sub;
  bool *owner_flag, *sub_flag;
  double *PhiB;
  double *weight, *ARinvCT, *CARinvCT, *lambda_r;

  int *amg_A1R, *amg_A2R, *amg_nodebegR, *amg_local_dofR, amg_nnodeR;
  double *amg_xR, *amg_yR, *amg_zR, *amg_matR;
  int *amg_A1I, *amg_A2I, *amg_nodebegI, *amg_local_dofI, amg_nnodeI;
  double *amg_xI, *amg_yI, *amg_zI, *amg_matI;
  int *amg_A1C, *amg_A2C, *amg_nodebegC, *amg_local_dofC, amg_nnodeC;
  double *amg_xC, *amg_yC, *amg_zC, *amg_matC;

  // old stuff
  solver_crd *AR, *AI, *AKc;
  int *mycdof, nmycon, cg_iter, n_orthog_used, ncon_global, max_added_corner;
  int nI, nB, nC, nR, nB_own, nBB, *dofI, *dofB, *dofC, *dofR, *dofBB, sub_begin;
  int num_tied_down, *tied_down;
  double *RHS_cg, *SOL_cg, *TEMP_cg, *SOL_Kc, *TEMP_Kc;
  double *rcurra, *rhoa, *betaa, *pApa, *Dtri, *Etri, *econa;
  double *P_ortho, *AP_ortho, *orth1, *orth2, *PhirTPhir;
  double *VV, *RR, *HH, *zz, *cc, *ss, *norms, *gmres_vec, *gmres_sum;
  ofstream fout;
};
#endif // CLIP_SOLVER2_HPP
