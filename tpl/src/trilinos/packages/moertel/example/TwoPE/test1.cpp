/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/

#include "Moertel_config.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include <assert.h>
#include "Epetra_MpiComm.h"
#include "mrtr_manager.H"
#include "mrtr_segment_linear1D.H"
#endif

//
// two dimensional example (1D interface) run on two processors
//

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  int i, MyPID, NumProc, nedge, side=0, nnode, ScanSum;
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
  NumProc = Comm.NumProc();
  if (NumProc != 2)
  {
    cout << "test1.exe needs to be run on 2 processors\n";
	return 1;
  }
  // ------------------------------------------------------------- //
  // create an empty MOERTEL::Interface, in this example just one
  // ------------------------------------------------------------- //
  int printlevel = 8; // ( moertel takes values 0 - 10 )
  MOERTEL::Interface interface(0, true, Comm, printlevel);
  // ------------------------------------------------------------- //
  // Add nodes on both sides of interface to interface
  // ------------------------------------------------------------- //
  if (MyPID == 0) {
    nedge = 4;
    side = 0;
  }
  if (MyPID == 1) {
    nedge = 5;
    side = 1;
  }
  nnode = nedge + 1;
  Comm.ScanSum(&nnode, &ScanSum, 1);
  int *nodeid = new int[nnode];
  double L = 1/double(nedge);
  double coord[3];
  coord[0] = 1;
  coord[2] = 0; 
  //
  // nodes on interface
  //
  for (i=0; i<nnode; i++) {
    coord[1] = i*L;
    nodeid[i] = ScanSum - nnode + i;
    MOERTEL::Node node(nodeid[i], coord, 1, &nodeid[i], false, printlevel);
    if(!interface.AddNode(node, side)){
    	cout << "interface AddNode returned false\n";
    	return 1;
	}
  }
  //
  // segments on interface
  //
  Comm.ScanSum(&nedge, &ScanSum, 1);
  int econ[2], seg_id;
  for (i=0; i<nedge; i++) {
    // segments have to be defined in anti-clockwise order for each subdomain!
    if (MyPID==0)
    {
      econ[0] = nodeid[i];
      econ[1] = nodeid[i+1];
    }
    else
    {
      econ[1] = nodeid[i];
      econ[0] = nodeid[i+1];
    }
    seg_id = ScanSum - nedge + i + 1;
    MOERTEL::Segment_Linear1D segment(seg_id, 2, econ, printlevel);
    if(!interface.AddSegment(segment, side)){
    	cout << "interface AddSegment returned false\n";
		return 1;
	}
  }
  //
  // mortar (master) side is 0
  //
  interface.SetMortarSide(0);
  interface.SetFunctionTypes(MOERTEL::Function::func_Linear1D,       // primal trace space
			     MOERTEL::Function::func_DualLinear1D);  // dual mortar space (recommended)
  if (!interface.Complete()) {
    cout << "Interface completion returned false\n";
	return -1;
  }
  // ------------------------------------------------------------- //
  // create an empty MOERTEL::Manager for 2D problems
  // It organizes everything from integration to solution
  // ------------------------------------------------------------- //
  MOERTEL::Manager manager(Comm, printlevel);
  manager.SetDimension(MOERTEL::Manager::manager_2D);
  // ------------------------------------------------------------- //
  // Add the interface to the manager
  // ------------------------------------------------------------- //
  manager.AddInterface(interface);
  
  // ------------------------------------------------------------- //
  // we have to supply a rowmap of the whole problem so Moertel knows
  // how the system of equations looks like. Constraint equations will 
  // have rows n to n + nlm-1 where n is the number of degrees of freedom
  // and nlm is the number of lagrange multipliers
  // ------------------------------------------------------------- //
  int gnnode;
  Comm.SumAll(&nnode, &gnnode, 1);
  Epetra_Map map(gnnode,nnode,nodeid,0,Comm);
  manager.SetProblemMap(&map);
  
  // ============================================================= //
  // choose integration parameters
  // ============================================================= //
  //Teuchos::ParameterList& moertelparams = manager.Default_Parameters();
  // this takes effect in 3D only
  //moertelparams.set("exact values at gauss points",true);
  // 1D interface possible values are 1,2,3,4,5,6,7,8,10 (2 recommended with linear shape functions)
  //moertelparams.set("number gaussian points 1D",2);
  // 2D interface possible values are 3,6,12,13,16,19,27 (12 recommended with linear functions)
  //moertelparams.set("number gaussian points 2D",12);

  // ============================================================= //
  // Here we are done with the construction phase of the interface
  // so we can integrate the mortar integrals
  // (Note we have not yet evaluated the PDE at all!)
  // ============================================================= //
  manager.Mortar_Integrate();
    
  // print interface information
  // (Manager, Interface, Segment, Node implement the << operator)
  if (printlevel) cout << manager;
  
  
  // get the 2 pieces of the constraint equation
  const Epetra_CrsMatrix* D = manager.D();
  const Epetra_CrsMatrix* M = manager.M();

  // create an empty Epetra_CrsMatrix and add D and M
  Epetra_CrsMatrix constraints(Copy,D->RowMap(),5);
  MOERTEL::MatrixMatrixAdd(*D,false,1.0,constraints,0.0);
  MOERTEL::MatrixMatrixAdd(*M,false,1.0,constraints,1.0);
  constraints.FillComplete(D->DomainMap(),D->RangeMap());
  constraints.OptimizeStorage();
  cout << constraints;
        
  MPI_Finalize();

  cout << "\nTest passed!" << endl;

  return 0;

#else

  return -1;

#endif

}
