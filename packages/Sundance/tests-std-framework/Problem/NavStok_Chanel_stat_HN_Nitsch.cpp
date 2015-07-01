/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */


#include "Sundance.hpp"
void CopyDiscreteFunction(Expr &copytothis, Expr &copythis)
{
// getvector gives a reference, setvector makes a shallow copy
  Vector<double> veccopythis = DiscreteFunction::discFunc(copythis)->getVector() ;
  Vector<double> veccopytothis = DiscreteFunction::discFunc(copytothis)->getVector() ;
  veccopytothis.acceptCopyOf(veccopythis) ;
  DiscreteFunction::discFunc(copytothis)->setVector(veccopytothis) ;
}

void CopyDiscreteFunctionaxpby(Expr &copyresulttothis, double a, Expr &x, double b, Expr &y)
{
// getvector gives a reference, setvector makes a shallow copy
  Vector<double> vecx = DiscreteFunction::discFunc(x)->getVector() ;
  Vector<double> vecy = DiscreteFunction::discFunc(y)->getVector() ;
  Vector<double> veccopyresulttothis = DiscreteFunction::discFunc(copyresulttothis)->getVector() ;
  veccopyresulttothis.acceptCopyOf(vecx);
  veccopyresulttothis=a*veccopyresulttothis+b*vecy;
  DiscreteFunction::discFunc(copyresulttothis)->setVector(veccopyresulttothis) ;
}


/** 
 * Solves the Navier-Stokes equation in 2D using Taylor-Hood-Elements:
 *
 * du/dt - nu* Delta u + (u*grad) u + grad p = f     in Omega x I
 *                                     div u = 0     in Omega x I
 *  Driven Cavity scenario
 */

/*
This code produces following results (with RefLevel = 1 and 220X41)
Drag: c_D(u,p) = 5.55722
Lift: c_L(u,p) = 0.0107113
*/

#define RefLevel 1

// Will be used to select the left, bottom, right, and top boundary
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.2) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-0.41) < 1.0e-10;})

REFINE_MESH_ESTIMATE(MeshRefEst , { \
    // this is the refinement function
        if ((cellPos[0] < 0.26) && (cellPos[0] > 0.14) && (cellPos[1] < 0.26) && (cellPos[1] > 0.14) && (cellLevel < RefLevel)) \
//           ||((cellPos[0] < 0.32) && (cellPos[0] > 0.12) && (cellPos[1] < 0.28) && (cellPos[1] > 0.12) && (cellLevel < 2)))
                                      return 1;\
                                 else return 0; } ,
  // this is the coarse load estimator
  { if (cellPos[0] < 0.7) {return 20;} else {return 1;} } )
  

MESH_DOMAIN( MeshDomain , {return true;})

int main(int argc, char** argv)
{
  
  try
	{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      int myrank = MPIComm::world().getRank();
 
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();
      
      /* refinement criterion object */
      RefinementClass refCl1 = new MeshRefEst();

      ParametrizedCurve curve_c = new Circle( 0.2 , 0.2 , 0.05 , 1 , 1e-7 );
      //ParametrizedCurve curve = new Polygon2D( "wing.txt" , 1.0 , 1e-7 );
            
      /* mesh domain estimation */
      MeshDomainDef meshDom_circl = new CurveDomain( curve_c , Outside_Curve );

      /* mesh domain, dummy class */
      MeshDomainDef meshDom = new MeshDomain();

      MeshType meshType = new HNMeshType2D();
      MeshSource mesher = new HNMesher2D(0.0, 0.0, 2.2 , 0.41 , 51 , 11 , meshType , refCl1 , meshDom_circl );
      Mesh mesh = mesher.getMesh();
      //curve.setMesh( mesh );
      
      // write circle to VTK 
      PolygonQuadrature::setNrMaxLinePerCell(6);
      ParametrizedCurve curve = curve_c.getPolygon(mesh, 0.07);
      curve.writeToVTK("p_wing.vtk");
            
      /* Object needed for curve integral calculation */
      ParametrizedCurve curveIntegral = new ParamCurveIntegral(curve);

      cout << "Nr Points  "<<mesh.numCells(0) << endl;
      cout << "My Rank is :" << myrank << endl;
      
      // Create cell filters
      CellFilter interior = new MaximalCellFilter();
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter edges = new DimensionalCellFilter(1);

      // Find left, right, top, bottom and cylinder boundary edges
      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
     
      /* the cell predicate for different integrations on curve */
      CellPredicate curveIN = new CellCurvePredicate( curve , Inside_Curve );
      CellPredicate curveOUT = new CellCurvePredicate( curve , Outside_Curve );
      CellPredicate curveON = new CellCurvePredicate( curve , On_Curve );
     
      CellFilter InCurve = interior.subset(curveIN);
      CellFilter OutsideCurve = interior.subset(curveOUT);
      CellFilter OnCurve = interior.subset(curveON);

      CellFilter curveFilt = edges.subset(curveIN);
      
       // Create unknown and test functions, discretized using Taylor Hood elements 
      BasisFamily L1 = new Lagrange(1);
      BasisFamily L2 = new Lagrange(2);
      Expr ux = new UnknownFunction(L2, "u_x");
      Expr vx = new TestFunction(L2, "v_x");
      Expr uy = new UnknownFunction(L2, "u_y");
      Expr vy = new TestFunction(L2, "v_y");
      Expr p = new UnknownFunction(L1, "p");
      Expr q = new TestFunction(L1, "q");
      Expr u = List(ux, uy);
      Expr v = List(vx, vy);

      // Create differential operator and coordinate functions 
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      // We need a quadrature rule for doing the integrations 
      QuadratureFamily quad2 = new GaussianQuadrature(4);
      QuadratureFamily quad4 = new GaussianQuadrature(6);
      QuadratureFamily quad_hi = new GaussLobattoQuadrature(6);
      QuadratureFamily quad_c12 = new GaussianQuadrature(12);
      QuadratureFamily quad_c = new PolygonQuadrature( quad_c12 );
      
      DiscreteSpace VelSpace(mesh, List(L2,L2), vecType);
      DiscreteSpace VelPreSpace(mesh, List(L2,L2,L1), vecType);

      // iterate from last timestep u_old 
      Expr uo = new DiscreteFunction(VelSpace, 0.0, "uo");
      // actual iterate (u,p) 
      Expr up = new DiscreteFunction(VelPreSpace, 0.0, "up");

      // zero right hand side 
      Expr fx = 0.0;
      Expr fy = 0.0;
      // viscosity parameter 
      Expr nu = new Sundance::Parameter(1.0/1000.0);
      Expr h = new CellDiameterExpr();
      Expr nx = new CurveNormExpr(0);
      Expr ny = new CurveNormExpr(1);
      Expr c1 = new UserDefOp(List(x,y), rcp(new CurveExpr(curve,0)) );
      curve.addNewScalarField( "fx" , 1.0 );
      curve.addNewScalarField( "fy" , 1.0 );
//      Expr nx = -(x-0.2)/sqrt((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2));
//      Expr ny = -(y-0.2)/sqrt((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2));
      double ux_Dirichlet = 0.0 , velx = 0.0; //0.04;
      double uy_Dirichlet = 0.0 , vely = 0.0; //0.01;
      double ga1=1000;
      double ga2=1000;
      double beta = 1e-5;

      Expr chix = 1.0/((1.0+exp(-200*(x-0.075)))*(1.0+exp(-15*(1.225-x))));
      Expr chiy = 1.0/((1.0+exp(-200*(y-0.075)))*(1.0+exp(-200*(0.33-y))));
      Expr chi=chix*chiy;

      // initial value 
      L2Projector projectorupar(VelSpace, List(0.0,0.0)); 
    
      L2Projector projectoruparp(VelPreSpace, List(0.0,0.0,0.0)); 
     
      // Define the Dirichlet BC 
      double hflow = 0.41;
      double Um = 0.3;
      double Umean = 2.0*(4.0*Um*((hflow-0.5*0.41)*0.5*0.41)/(hflow*hflow))/3.0;
      double rho = 1.0;
      double D = 0.1;
      double mu = 0.001;
      
      Expr uInflow = 4*Um*((hflow-y)*y)/(hflow*hflow);
      // if we add the edges then we can make sure that no fluid is escaping
      CellFilter walls = bottom + top; //+ curveFilt;
      Expr bc = EssentialBC(walls , (1.0/h)*v*u, quad4 )  //curveFilt
              + EssentialBC(left, (1.0/h)*(vx*(ux-uInflow)+vy*uy), quad4 );

/*
      Expr un = ux*nx+uy*ny;
      Expr unm = 0.5*(un-fabs(un));
      Expr udln = ux_Dirichlet*nx+uy_Dirichlet*ny;
      Expr udlnm = 0.5*(udln-fabs(udln));
*/
      Expr unm=1e-12*x;
      Expr udlnm=1e-12*x;
      Expr eqn = Integral(OutsideCurve, nu*(grad*vx)*(grad*ux)  
                          + nu*(grad*vy)*(grad*uy)
                          + vx*(ux*(dx*ux)+uy*(dy*ux))
                          + vy*(ux*(dx*uy)+uy*(dy*uy))
			  + 1e-2 * (grad*p)*(grad*q)
                          // symmetrization
//                          +0.5*(dx*ux+dy*uy)*(vx*ux+vy*uy)
                          - p*(dx*vx+dy*vy)
                          + q*(dx*ux+dy*uy),
                          quad4)
                 + Integral(OnCurve, nu*(grad*vx)*(grad*ux)  
                          + nu*(grad*vy)*(grad*uy)
                          + vx*(ux*(dx*ux)+uy*(dy*ux))
                          + vy*(ux*(dx*uy)+uy*(dy*uy))
			  + 1e-2 * (grad*p)*(grad*q)
                          // symmetrization
//                          +0.5*(dx*ux+dy*uy)*(vx*ux+vy*uy)
                          - p*(dx*vx+dy*vy)
                          + q*(dx*ux+dy*uy),
                          quad_hi , curve)
                 + Integral(OnCurve,
                            -nu*(dx*ux*nx+dy*ux*ny)*vx
                            -nu*(dx*uy*nx+dy*uy*ny)*vy
                            -nu*(dx*vx*nx+dy*vx*ny)*(ux-ux_Dirichlet)
                            -nu*(dx*vy*nx+dy*vy*ny)*(uy-uy_Dirichlet)
                            +p*(nx*vx+ny*vy)
                            -q*(nx*(ux-ux_Dirichlet)+ny*(uy-uy_Dirichlet))
                            +nu*ga1/h*((ux-ux_Dirichlet)*vx+(uy-uy_Dirichlet)*vy)
                            +ga2/h*((ux-ux_Dirichlet)*nx+(uy-uy_Dirichlet)*ny)
                            *(vx*nx+vy*ny)
                            -unm*(ux*vx+uy*vy)
                            +udlnm*(ux_Dirichlet*vx+uy_Dirichlet*vy),
                            quad_c , curveIntegral )
// New
// good: 1e-7, also OK: 1e-6
                   + Integral(InCurve+OnCurve, 1e-7*nu*(grad*vx)*(grad*ux)+1e-7*vx*ux
                                      +1e-7*nu*(grad*vy)*(grad*uy)+1e-7*vy*uy+1e-7*p*q,quad4)
                   - Integral(OnCurve, 1e-7*nu*(grad*vx)*(grad*ux)+1e-7*vx*ux
                                      +1e-7*nu*(grad*vy)*(grad*uy)+1e-7*vy*uy+1e-7*p*q,quad_hi, curve);
      // Create a TSF NonlinearOperator object 
      NonlinearProblem prob(mesh, eqn , bc, List(v, q),
                               List(u, p), up, vecType );

      ParameterXMLFileReader reader("nox-amesos.xml"); //nox-NavStok.xml
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);

      // Solve the nonlinear system
      SolverState<double> status = prob.solve(solver);

      DiscreteSpace discreteSpace(mesh, List( L1 , L1 ), vecType);
      L2Projector projector(discreteSpace, grad*up[2]);
      Expr p_grad = projector.project();
      
      DiscreteSpace discreteSpaceUX(mesh, List( L2 , L2 ), vecType);
      L2Projector projectorUX(discreteSpaceUX, grad*up[0]);
      Expr ux_grad = projectorUX.project();
      
      DiscreteSpace discreteSpaceUY(mesh, List( L2 , L2 ), vecType);
      L2Projector projectorUY(discreteSpaceUY, grad*up[1]);
      Expr uy_grad = projectorUY.project();
            
      // Write the field in VTK format
      FieldWriter w = new VTKWriter("NavStok_Benchmark");
      w.addMesh(mesh);

      // this is what should be passed for vector field plotting
      Expr expr_vector(List(up[0],up[1]));
      Expr ux_grad_E(List(ux_grad[0],ux_grad[1]));
      Expr uy_grad_E(List(uy_grad[0],uy_grad[1]));
      Expr expr_p_grad(List(p_grad[0],p_grad[1]));

      w.addField("ux", new ExprFieldWrapper(up[0]));
      w.addField("ux_grad", new ExprFieldWrapper(ux_grad_E));
      w.addField("uy", new ExprFieldWrapper(up[1]));
      w.addField("uy_grad", new ExprFieldWrapper(uy_grad_E));
      w.addField("p_grad", new ExprFieldWrapper(expr_p_grad));
      w.addField("p_grad_x", new ExprFieldWrapper(p_grad[0]));
      w.addField("p_grad_y", new ExprFieldWrapper(p_grad[1]));

      // add a vector field , try pnly one at once
      w.addField("vel", new ExprFieldWrapper(expr_vector));
      w.addField("p", new ExprFieldWrapper(up[2]));
      w.write();

      double tol = 1.0e-3;
      Expr errExpr1 = Integral(right , (up[0]-uInflow)*(up[0]-uInflow) , quad4);
      FunctionalEvaluator errInt1(mesh, errExpr1);
      double errorSq = errInt1.evaluate();
      double err_i = std::sqrt(errorSq);
      Out::os() << " error ux =" << err_i << endl;
      Sundance::passFailTest( err_i , tol);
      
      Expr errX = up[0]-ux_Dirichlet;
      Expr errExprX = Integral( OnCurve , errX*errX , quad_c , curveIntegral );
      FunctionalEvaluator errIntX( mesh , errExprX );
      double errorSqX = errIntX.evaluate();
      Out::os() << "Curve Integration X error = " << errorSqX << "  , sqrt(error)=" << sqrt(errorSqX) <<std::endl;
      Sundance::passFailTest( sqrt(errorSqX) , 0.01);
      
      Expr errY = up[1]-uy_Dirichlet;
      Expr errExprY = Integral( OnCurve , errY*errY , quad_c , curveIntegral );
      FunctionalEvaluator errIntY( mesh , errExprY);
      double errorSqY = errIntY.evaluate();
      Out::os() << "Curve Integration Y error = " << errorSqY << "  , sqrt(error)=" << sqrt(errorSqY) <<std::endl;
      Sundance::passFailTest( sqrt(errorSqY) , 0.01);
      
      Expr tmpExpr = Integral( OnCurve ,up[2] , quad_c , curveIntegral ); //, watchMe );
      FunctionalEvaluator tmpInt( mesh , tmpExpr );
      //double tmpVal = tmpInt.evaluate();
      //Out::os() << "tmpVal = " << tmpVal << std::endl;

      Expr dragExpr_V = Integral( OnCurve ,
	  - rho*mu*ny*( nx*dx*(ny*up[0] - nx*up[1]) + ny*dy*(ny*up[0] - nx*up[1]) ) , quad_c , curveIntegral );
      FunctionalEvaluator dragInt_V( mesh , dragExpr_V );
      double dragVal_Vel = dragInt_V.evaluate();
      Expr dragExpr_P = Integral( OnCurve , nx*up[2] , quad_c , curveIntegral );
      FunctionalEvaluator dragInt_P( mesh , dragExpr_P );
      double dragVal_Pres = dragInt_P.evaluate();
      double dragVal = dragVal_Vel + dragVal_Pres;
      Out::os() << "Drag_V = " << 2*dragVal_Vel/(rho*Umean*Umean*D) << std::endl;
      Out::os() << "Drag_P = " << 2*dragVal_Pres/(rho*Umean*Umean*D) << std::endl;
      Out::os() << "Drag = " << 2*dragVal/(rho*Umean*Umean*D) << std::endl;

      // =========== POLYGON EVAL ==============================
      Polygon2D* polyg = dynamic_cast<Polygon2D*>(curve.ptr().get());
      ParametrizedCurve curve_twin = polyg->createTwinPolygon( 0.1 , 0.1 , 0.5 , 0.5 );
      curve_twin.setMesh(mesh);
      
      Expr dragExpr_Curve = Integral( OnCurve , - rho*mu*ny*( nx*dx*(ny*up[0] - nx*up[1]) + ny*dy*(ny*up[0] - nx*up[1]) )
                    + nx*up[2] , quad_c , curveIntegral );
      FunctionalEvaluator dragInt_Curve( mesh , dragExpr_Curve );
      Expr liftExpr_Curve = Integral( OnCurve , rho*mu*nx*( nx*dx*(ny*up[0] - nx*up[1]) + ny*dy*(ny*up[0] - nx*up[1]))
                    + ny*up[2], quad_c , curveIntegral );
      FunctionalEvaluator liftInt_Curve( mesh , liftExpr_Curve );
      curve.setSpaceValues( dragInt_Curve , 0 );
      curve.setSpaceValues( liftInt_Curve , 1 );
      curve.writeToVTK("p_wing_withvalue.vtk");
      curve_twin.writeToVTK("p_wing_twin_withvalue.vtk");
      // =========== END POLYGON EVAL ==============================

      Expr liftExpr_V = Integral( OnCurve ,
	  rho*mu*nx*( nx*dx*(ny*up[0] - nx*up[1]) + ny*dy*(ny*up[0] - nx*up[1])) , quad_c , curveIntegral ); 
      FunctionalEvaluator liftInt_V( mesh , liftExpr_V );
      double liftVal_V = liftInt_V.evaluate();
      Expr liftExpr_P = Integral( OnCurve , ny*up[2] , quad_c , curveIntegral ); 
      FunctionalEvaluator liftInt_P( mesh , liftExpr_P );
      double liftVal_P = liftInt_P.evaluate();
      double liftVal = liftVal_V + liftVal_P;
      Out::os() << "Lift_V = " << 2*liftVal_V/(rho*Umean*Umean*D) << std::endl;
      Out::os() << "Lift_P = " << 2*liftVal_P/(rho*Umean*Umean*D) << std::endl;
      Out::os() << "Lift = " << 2*liftVal/(rho*Umean*Umean*D) << std::endl;
      
// -------------- Other lift and drag coefficient calculation ----------------
      Expr ncx=-nx;
      Expr ncy=-ny;
      Expr tcx=ncy;
      Expr tcy=-ncx;
      Expr FDragExpr = Integral(OnCurve,
                          nu*((dx*up[0]*tcx+dx*up[1]*tcy)*ncx
                              +(dy*up[0]*tcx+dy*up[1]*tcy)*ncy)*ncy
                          -up[2]*ncx,
                              quad_c , curveIntegral);
      double c_D = 2.0*evaluateIntegral(mesh,FDragExpr)/(Umean*Umean*D);
      cerr << "Drag: c_D(u,p) = " << c_D << endl << endl;
      Sundance::passFailTest( ::fabs(c_D - 5.57) , 1.0);

      Expr Phix=chi;
      Expr Phiy=0.0;

      Expr FDragExprOm = Integral(OutsideCurve,
                          - Phix*(up[0]*(dx*up[0])+up[1]*(dy*up[0]))
                          - Phiy*(up[0]*(dx*up[1])+up[1]*(dy*up[1]))
                          // symmetrization
////                          -0.5*(dx*up[0]+dy*up[1])*(Phix*up[0]+Phiy*up[1])
                          - nu*(2.0*(dx*up[0])*(dx*Phix)
                               +(dy*up[0]+dx*up[1])*(dy*Phix+dx*Phiy)
                               +2.0*(dy*up[1])*(dy*Phiy))
                          + up[2]*(dx*Phix+dy*Phiy),
                          quad4)
                 + Integral(OnCurve,
                          - Phix*(up[0]*(dx*up[0])+up[1]*(dy*up[0]))
                          - Phiy*(up[0]*(dx*up[1])+up[1]*(dy*up[1]))
                          // symmetrization
////                          -0.5*(dx*up[0]+dy*up[1])*(Phix*up[0]+Phiy*up[1])
                          - nu*(2.0*(dx*up[0])*(dx*Phix)
                               +(dy*up[0]+dx*up[1])*(dy*Phix+dx*Phiy)
                               +2.0*(dy*up[1])*(dy*Phiy))
                          + up[2]*(dx*Phix+dy*Phiy),
                          quad_hi , curve);
      double c_D_Om = 2.0*evaluateIntegral(mesh,FDragExprOm)/(Umean*Umean*D);
      cerr << "Drag: c_D_Om(u,p) = " << c_D_Om << endl << endl;
      Sundance::passFailTest( ::fabs(c_D_Om - 5.57) , 1.0);

      Expr FLiftExpr = Integral(OnCurve,
                          -nu*((dx*up[0]*tcx+dx*up[1]*tcy)*ncx
                               +(dy*up[0]*tcx+dy*up[1]*tcy)*ncy)*ncx
                          -up[2]*ncy,
                              quad_c , curveIntegral);

      double c_L = 2.0*evaluateIntegral(mesh,FLiftExpr)/(Umean*Umean*D);

      cerr << "Lift: c_L(u,p) = " << c_L << endl << endl;
      Sundance::passFailTest( ::fabs(c_L - 0.0107) , 0.2);

      Expr Psix=0.0;
      Expr Psiy=chi;

      Expr FLiftExprOm = Integral(OutsideCurve,
                          - Psix*(up[0]*(dx*up[0])+up[1]*(dy*up[0]))
                          - Psiy*(up[0]*(dx*up[1])+up[1]*(dy*up[1]))
                          // symmetrization
//                          -0.5*(dx*up[0]+dy*up[1])*(Psix*up[0]+Psiy*up[1])
                          - nu*(2.0*(dx*up[0])*(dx*Psix)
                               +(dy*up[0]+dx*up[1])*(dy*Psix+dx*Psiy)
                               +2.0*(dy*up[1])*(dy*Psiy))
                          + up[2]*(dx*Psix+dy*Psiy),
                          quad4)
                 + Integral(OnCurve,
                          - Psix*(up[0]*(dx*up[0])+up[1]*(dy*up[0]))
                          - Psiy*(up[0]*(dx*up[1])+up[1]*(dy*up[1]))
                          // symmetrization
//                          -0.5*(dx*up[0]+dy*up[1])*(Psix*up[0]+Psiy*up[1])
                          - nu*(2.0*(dx*up[0])*(dx*Psix)
                               +(dy*up[0]+dx*up[1])*(dy*Psix+dx*Psiy)
                               +2.0*(dy*up[1])*(dy*Psiy))
                          + up[2]*(dx*Psix+dy*Psiy),
                          quad_hi , curve);
      double c_L_Om = 2.0*evaluateIntegral(mesh,FLiftExprOm)/(Umean*Umean*D);
      cerr << "Lift: c_L_Om(u,p) = " << c_L_Om << endl << endl;
      Sundance::passFailTest( ::fabs(c_L_Om - 0.0107) , 0.2);
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
  return(0);
}
