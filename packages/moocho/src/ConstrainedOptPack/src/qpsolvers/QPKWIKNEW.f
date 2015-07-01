*********************************
*		QPKWIKNEW Interface		*
*********************************

*=============================================================================================

      SUBROUTINE QPKWIKNEW (
*                *** Input
     &              N,M1,M2,M3,GRAD,UINV_AUG,LDUINV_AUG,IBND,BL,BU
     &                  ,A,LDA,YPY,IYPY,WARM,NUMPARAM,MAX_ITER
*                *** Input/Output
     &              ,X,NACTSTORE,IACTSTORE,INF
*                *** Output
     &              ,NACT,IACT,UR,EXTRA,ITER,NUM_ADDS,NUM_DROPS
*                *** Internal state
     &              ,ISTATE
*                *** Workspace
     &              ,LRW,RW     )

* Need documentation!
*
* INF out = ( -1: constraints inconsistant, -2: LRW too small )
*
*
* QPKWIKNEW  Solves a quadratic programming (QP) subproblem in an rSQP algorithm or a
* isolated QP.  For this solver, the second order matrix must be positive definate
* so that the QP is convex.  This is because the second order matrix is passed in
* as the inverse of the cholesky factor and this factor only exists if this
* (symmetric) matrix is positive definate.
* This QP solver was designed primarily for use in solving a reduced
* rSQP subproblem and a relaxation is designed with that in mind.  With that
* being the case, the description of QPKWIKNEW will be given in two seperate sections:
* one for solving general QP's and one for solving an rSQP reduced QP subprolem.  This
* will result in more documentation but should also be more understandable.  Each section
* has its own mathematical nomenclature and mapping to the arguments for QPKWIKNEW.
*
* ----------------------------------------
* (1) QPKWIKNEW for solving a general QP
* ----------------------------------------
* 
* The form of this (relaxed) QP is given below (In MATLAB notation):
*
*	min		g'*x + (1/2)*x'*H*x + (eta + eta^2/2)*M
*
*	s.t.	xL			<=	x						<= xU
*			eL			<=	E*x						<= eU
*							F*x + eta*f				 = f
*			0			<=	eta						<= 1
*
* The Lagrangian for this QP is:
*
*	L = g'*x + (1/2)*x'*H*x + (eta + eta^2/2)*M
*       + nuL' * (xL - x)
*	    + nuU' * (x - xU)
*       + gamaL' * (eL - E*x)
*	    + gamaU' * (E*x + eU)
*       + lambda' * (F*x - (1-eta)*f)
*		+ kapaL' * (0 - eta)
*	    + kapaU' * (eta - 1)
*
* The optimality conditions (first order KKT) for this QP are:
*
*	Linear dependance of gradients:
*	-------------------------------
*	del(L,x)   = g + H*x + nu + E'*gama + F'*lambda = 0
*	d(L)/d(eta) = (1+eta) * M + lambda'*f + kapa = 0
*	            where :
*					nu = nuU - nuL
*					gama = gamaU - gamaL
*					kapa = kapaU - kapaL
*
*	Feasibility:
*	-----------
*	xL			<=	x						<= xU
*	eL			<=	E*x						<= eU
*					F*x + eta*f				 = f
*	0			<=	eta						<= 1
*
*	Complementarity:
*	---------------
*	nuL' * (xL - x) = 0
*	nuU' * (x - xU) = 0
*	gamaL' * (eL - E*x) = 0
*	gamaU' * (E*x + eU) = 0
*	kapaL' * (0 - eta) = 0
*	kapaU' * (eta - 1) = 0
*
*	Multiplier non-negativity:
*	-------------------------
*	nuL, nuU, gamaL, gamaU, kapaL, kapaU >= 0
*
*	In this context, the relaxation parameter eta represents a relaxation of the equality
*	constraints.  If the original QP has no solution (eta > 0) then the equality constraints
*	are relaxed until a feasible region is found.  The parameter M is made large enough
*	so that if a feasible region exists then the solution will be found within it, otherwise
*	the feasible region will only be opened as little as needed.
*
*	The following is the mapping of the above QP into the parameter names for QPKWIKNEW
*	(in MATLAB notation): 
*
*	min		GRAD'*X + (1/2)*X'*inv(UINV*UINV')*X + (EXTRA + EXTRA^2 / 2) * M
*
*	s.t.	BL(1:M1)		<= X(i)									<= BU(1:M1)
*			BL(M1+1:M1D)	<= A(1:M2,:)*X							<= BU(M1+1:M1D)
*							   A(M2+1:M2D)*X						 = BU(M2D+1:M3D)
*			0.0				<= EXTRA								<= 1.0
*
*	 Where:
*		i is in the set IBND(1:M1) corresponding to variables with bounds
*
*---------------------------------------------------------------------------------------------
*		Formal Parameter Descriptions
*---------------------------------------------------------------------------------------------
*
*	Name		I/O		Type	Meaning
*
*	N			I		I		Dimension of X
*	M1			I		I		Number of variables with simple bounds
*	M2			I		I		Number of rows in E
*	M3			I		I		Number of rows in F
*	GRAD(N)		I		R*8		First order objective vector g
*	UINV_AUG(	I		R*8		Inverse of the upper Cholesky factor of second order objective
*		LDUINV_AUG,N+1)			matrix H augmented with the terms for the relaxation parameter.
*								See the subroutine QPKWIKNEW_INIT_UINV_AUG(...).
*	LDUINV_AUG	I		I		The leading dimension of UINV_AUG
*	IBND(M1)	I		I		Contains the indices of the M1 bounded variables in X:
*									IBND(1:M1) = { i | xL(i) > -inf OR xU(i) < inf }
*	BL(M1+M2)	I		R*8		The first M1 elements are the lower bounds on variables X:
*									BL(i) = xL(IBND(i)), for i = 1,...,M1
*								The last M2 elements are the lower bounds on
*								the inequality constraints:
*									BL(i+M1) = eL(i), for i = 1,...,M2
*	BU(			I		R*8		The first M1 elements are the upper bounds on the bounded
*	M1_M2_M3)					variables x:
*									BU(i) = xU(IBND(i)), for i = 1,...,M1
*								The next M2 elements are the upper bounds on the
*								inequality constraints:
*									BU(i+M1) = eU(i), for i = 1,...,M2
*								The last M3 elements are the right hand sides of the
*								equality constriants:
*									BU(M1+M2+1:M1+M2+M3) = f
*	A(LDA,N)	I		R*8		The first M2 rows is the coefficient matrix for the 
*								linear inequality constraints matrix:
*									A(1:M2,:) = E
*								The last M3 rows is the coeffient matrix for the linear
*								equality constraints matrix:
*									A(M2+1:M2+M3,:) = F
*	LDA			I		I		The leading dimension of A
*	YPY(M1+M2)	I		R*8		Set all elements to 0.0
*	IYPY		I		I		This is an undocumented parameter that as far as I can
*								tell is always set to 1 and if it is set to zero then the
*								sign of an input element in YPY is changed for some condition
*								that I don't understand and is not documented.  Thanks.
*	WARM		I		I		If WARM = 1 then QPKWIKNEW will use a warm start based on the
*								previous call.  The algorithm will add violated constraints
*								contained in IACTSTORE first.
*   NUMPARAM(3)	I		R*8		Numerical parameters used in the computations.  Their meaning
*								and default values are:
*									NUMPARAM(1)	= SMALL		= 1e-10
*									NUMPARAM(2)	= VSMALL	= 1e-20
*									NUMPARAM(3)	= VLARGE	= 1e+20  (Used as M)
*   MAX_ITER    I       I       Maxinum number of allowed QP iterations
*	X(N)		O		R*8		The solution vector x.
*	NACTSTORE	I/O		I		The number of active constraints used for a warm start
*								on input for a warm start and is set on output to be
*								used on the next iteration.
*	IACTSTORE(	I/O		I		On input for a warm start it gives the indices of the
*		N)						guess of the active constraints and on output gives
*								the active constraints of the solution.  This array
*								is set by this subroutine and need not be accessed
*								by the user.  The values in IACTSTORE(*) have the
*								same meaning as in IACT(*).
*	INF			I/O		I		Status/control flag.  Set to zero for the first call
*								and for any call where the sparsity pattern for A
*								has changed from the last call.
*								On output it gives the error condition.
*									>= 0	Success
*									-1		constraints inconsistant
*									-2		LRW too small
*                                   -3      MAX_ITER exceeded
*	NACT		O		I		Gives the number of active constraints at the solution.
*	IACT(N+1)	I/O		R*8		Gives the indices of active constraints (non-zero multipliers).
*								On output gives the active set.
*
*						[1			, M1		  ]	-> nuL(IBND(j))
*						[M1+1		, M1+M2		  ]	-> gamaL(j-M1)
*		j = IACT(i) =	[M1+M2+1	, 2*M1+M2	  ]	-> nuU(IBND(j-M1-M2))
*						[2*M1+M2+1	, 2*M1+2*M2	  ]	-> gamaU(j-2*M1-M2)
*						[2*M1+2*M2+1, 2*M1+2*M2+M3]	-> lambda(j-2*M1-2*M2)
*						[	 2*M1+2*M2+M3+1		  ]	-> +- kapa	
*
*	UR(N+1)		O		R*8		The non-zero Lagrange multpliers for the QP constraints.
*								UR(i) is the Lagrange multplier for the active constraint
*								IACT(i).  These multipliers are categorized according to
*								IACT(i) (see IACT).
*	EXTRA		O		R*8		The slack parameter eta [0,1].  EXTRA = 0 for a feasible QP.
*	ITER		O		I		Gives the number of QP iterations used.
*   NUM_ADDS    O       I       Gives the number of QP iterations were a constraint was added
*   NUM_DROPS   O       I       Gives the number of QP iterations were a constraint was dropped
*	ISTATE(*)	S		I		Array for internal state data.  The length this array
*								must be at least as large as the value returned from
*								the function QPKWIKNEW_LISTATE(...)
*	LRW			W		I		Length of real workspace.  This length is returned
*								by the function QPKWIKNEW_LRW(...).
*	RW(LRW)		W		R*8		Real workspace
*
* -----------------------------------------------------------------------
* (2) QPKWIKNEW for solving a reduced QP subproblem in an rSQP algorithm
* -----------------------------------------------------------------------
*
* The form of this reduced QP subproblem (relaxed) is given below (In MATLAB notation):
*
*	min		rGf' * pz + (1/2) * pz' * rHL * pz + (eta + eta^2/2) * M
*
*	s.t.
*	xL(indep) - xk(indep) - Ypy(indep)	<=	pz - eta*Ypy(indep)	<= xU(indep) - xk(indep) - Ypy(indep)
*	xL(dep) - xk(dep) - Ypy(dep)		<=	D*x - eta*Ypy(dep)	<= xU(dep)  - xk(dep) - Ypy(dep)
*										    V*x					 = -(U*py + c(dep))
*	0									<=	eta					<= 1
*
* See the rSQP++ nomenclatrue guide and C++ name mapping for a description of these quantities.
* The equality constriants are not relaxed as they technically should be.
* This QP is stated for a variable reduction decomposition for Z.
* The original QP inequality constraint (for the relaxed full QP
* , d = Z * pz + (1-eta) * Y * py) is:
*
*	xL	<=	xk + Z * pz + (1-eta)*Ypy	<= xU
*
* The variable reduction decomposition partitions x = [ x(dep)' ; x(indep)' ]' and Z is:
*
*	Z = [ D' ; I ]	where D = - inv(C) * N
*
* The vectors xL, xU, and Ypy also become partitioned into dependent and independent vectors.
* Applying this variable reduction Z yields the two inequality constriants shown.
*
* The Lagrangian for this QP is:
*
*	L = rGf'*pz + (1/2)*pz'*rHL* pz + (eta + eta^2/2)*M
*       + nuL(indep)' * (xL(indep) - xk(indep) - pz - (1-eta)*Ypy(indep))
*	    + nuU(indep)' * (pz + (1-eta)*Ypy(indep) - xU(indep) + xk(indep))
*       + nuL(dep)' * (xL(dep) - xk(dep) - D*pz - (1-eta)*Ypy(dep))
*	    + nuU(dep)' * (D*pz + (1-eta)*Ypy(dep) - xU(dep) + xk(dep))
*       + lambda(dep)' * (V*x + (U*py + c(dep)))
*		+ kapaL' * (0 - eta)
*	    + kapaU' * (eta - 1)
*
* The optimality conditions (first order KKT) for this QP are:
*
*	Linear dependance of gradients:
*	-------------------------------
*	del(L,pz)   = rGf + rHL*pz + nu(indep) + D'*nu(dep) + V'*lambda(dep) = 0
*
*	d(L)/d(eta) = (1+eta) * M - nu(indep)'*Ypy(indep) - nu(dep)'*Ypy(dep)
*	              + kapa = 0
*
*				where :
*					nuL = [ nuL(indep) ; nuL(dep) ]
*					nuU = [ nuU(indep) ; nuU(dep) ]
*					nu = nuU - nuL
*					kapa = kapaU - kapaL
*
*	Feasibility:
*	-----------
*	xL(indep) - xk(indep)	<=	pz + (1-eta)*Ypy(indep)		<= xU(indep) - xk(indep)
*	xL(dep) - xk(dep)		<=	D*x + (1-eta)*Ypy(dep)		<= xU(dep)  - xk(dep)
*								V*x							 = -(U*py + c(dep))
*	0						<=	eta							<= 1
*
*	Complementarity:
*	---------------
*	nuL(indep)' * (xL(indep) - xk(indep) - pz - (1-eta)*Ypy(indep))
*	nuU(indep)' * (pz + (1-eta)*Ypy(indep) - xU(indep) + xk(indep))
*	nuL(dep)' * (xL(dep) - xk(dep) - D*pz - (1-eta)*Ypy(dep))
*	nuU(dep)' * (D*pz + (1-eta)*Ypy(dep) - xU(dep) + xk(dep))
*	kapaL' * (0 - eta)
*	kapaU' * (eta - 1)
*
*	Multiplier non-negativity:
*	-------------------------
*	nuL, nuU, kapaL, kapaU >= 0
*
*	This QP is structured in such a way as to be useful in solving a reduced QP subproblem of an
*	rSQP algorithm. The variable eta is added to allow a relaxation of the feasible region and
*	allow a solution.  If no other solution is possible, this QP will always have the solution
*	pz = 0.0, eta = 1.0 if there are no equality constraints.  The value of the constant M is made
*   large enough so that if there is
*	a feasible solution then eta will be driven to zero.  If there is not a feasible solution
*	then eta will be greater than zero.  In this way the rSQP step is cut back in such a
*	way that the bounds on the variables will never be violated.  See Schmid and Biegler (1994)
*	"Quadratic Programming Methods for Reduced Hessian SQP" in "comp. in chem. engng" vol. 18
*	pp. 817 - 832 for a full descrition of the method. 
*
*	The following the mapping of the above reduced QP into the parameter names for QPKWIKNEW
*	(in MATLAB notation): 
*
*	min	GRAD' * x + (1/2) * X' * inv(UINV * UINV') * X + (EXTRA + EXTRA^2 / 2) * M
*
*	 s.t.   BL(1:M1)		<= X(i) - EXTRA * YPY(1:M1)					<= BU(1:M1)
*			BL(M1+1:M1D)	<= A(1:M2,:) * X - EXTRA * YPY(M1+1:M1+M2)	<= BU(M1+1:M1D)
*							   A(M2+1:M2D,:) * X						 = BU(M1D+1:M3D)
*			0.0				<= EXTRA									<= 1.0
*
*	 Where:
*		i is in the set IBND(1:M1) corresponding to variables with bounds
*
*---------------------------------------------------------------------------------------------
*		Formal Parameter Descriptions
*---------------------------------------------------------------------------------------------
*
*	Name		I/O		Type	Meaning
*
*	N			I		I		Dimension of pz
*	M1			I		I		Number of variables with simple bounds
*	M2			I		I		Number of rows in V
*	M3			I		I		Number of rows in U
*	GRAD(N)		I		R*8		First order objective vector rGf
*	UINV_AUG(	I		R*8		Inverse of the upper Cholesky factor of second order objective
*		LDUINV_AUG,N+1)			matrix rHL augmented with the terms for the relaxation parameter.
*								See the subroutine QPKWIKNEW_INIT_UINV_AUG(...).
*	LDUINV_AUG	I		I		The leading dimension of UINV_AUG
*	IBND(M1+M2)	I		I		The first M1 elements contains the indices of the bounded
*								variables in pz:
*									IBND(1:M1) = { i | xL(indep)(i) > -inf AND xU(indep)(i) < inf }
*								The last M2 elements are the indices of the bounded dependent full
*								space variables.  Strictly speaking, these elements are not
*								needed by QPKWIKNEW but are included because of the original Fortran
*								implemenation.  They are used to help identify Lagrange multipliers
*								latter:
*									IBND(M1+1:M1+M2) = { i | xL(dep)(i) > -inf AND xU(dep)(i) < inf }
*	BL(M1+M2)	I		R*8		The first M1 elements are the lower bounds on variables pz:
*									BL(i) = xL(IBND(i)) - xk(IBND(i)) - Ypy(IBND(i)), for i = 1,...,M1
*								The last M2 elements are the lower bounds on
*								the inequality constraints:
*									BL(i+M1) = xL(IBND(i+M1)) - xk(IBND(i+M1)) - Ypy(IBND(i)), for i = 1,...,M2
*	BU(			I		R*8		The first M1 elements are the upper bounds on the bounded
*		M1+M2+M3 )				variables pz:
*									BU(i) = xU(IBND(i)) - xk(IBND(i)) - Ypy(IBND(i)), for i = 1,...,M1
*								The next M2 elements are the upper bounds on the
*								inequality constraints:
*									BU(i+M1) = xU(IBND(i+M1)) - xk(IBND(i+M1)) - Ypy(IBND(i)), for i = 1,...,M2
*								The last M3 elements are the right hand sides of the
*								equality constriants:
*									BU(M1D+1:M3D) = -(U*py + c(dep))
*	A(LDA,N)	I		R*8		The first M2 rows is the coefficient matrix for the 
*								dense part of Z, D = - inv(C) * N:
*									A(1:M2,:) = D
*								The last M3 rows is the coeffient matrix for the
*								dependent constriant matrix V (see rSQP++ guide):
*									A(M2+1:M2+M3,:) = V
*								linear inequality constraints matrix:
*									A(1:M2,:) = E
*								The last M3 rows is the coeffient matrix for the linear
*								equality constraints matrix:
*									A(M2+1:M2+M3,:) = F
*	LDA			I		I		The leading dimension of A
*	YPY(M1+M2)	I		R*8		Elements of Ypy for the bounded variables,  Specifically:
*									YPY(i) = Ypy(IBND(i)), for i = 1,...,M1+M2.
*	IYPY		I		I		This is an undocumented parameter that as far as I can
*								tell is always set to 1 and if it is set to zero then the
*								sign of an input element in YPY is changed for some condition
*								that I don't understand and is not documented.  Thanks.
*	WARM		I		I		If WARM = 1 then QPKWIKNEW will use a warm start based on the
*								previous call.  The algorithm will add violated constraints
*								contained in IACTSTORE first.
*   NUMPARAM(3)	I		R*8		Numerical parameters used in the computations.  Their meaning
*								and default values are:
*									NUMPARAM(1)	= SMALL		= 1e-10
*									NUMPARAM(2)	= VSMALL	= 1e-20
*									NUMPARAM(3)	= VLARGE	= 1e+20  (Used as M)
*	X(N)		O		R*8		The solution vector pz.
*	NACTSTORE	I/O		I		The number of active constraints used for a warm start
*								on input for a warm start and is set on output to be
*								used on the next iteration.
*	IACTSTORE(	I/O		I		On input for a warm start it gives the indices of the
*		N)						guess of the active constraints and on output gives
*								the active constraints of the solution.  This array
*								is set by this subroutine and need not be accessed
*								by the user.  The values in IACTSTORE(*) have the
*								same meaning as in IACT(*).
*	INF			I/O		I		Status/control flag.  Set to zero for the first call
*								and for any call where the sparsity pattern for A
*								has changed from the last call.
*								On output it gives the error condition.
*									>= 0	Success
*									-1		constraints inconsistant
*									-2		LRW too small
*                                   -3      MAX_ITER exceeded
*	NACT		O		I		Gives the number of active constraints at the solution.
*	IACT(N+1)	I/O		R*8		Gives the indices of active constraints (non-zero multipliers).
*								On output gives the active set.
*
*	IACT(N+1)	I/O		R*8		Gives the indices of active constraints (non-zero multipliers).
*								If WARM == 1 then on input the active constraints
*								stored in IACT are used as the initial guess for the
*								active set.  On output gives the active set.
*
*						[1			, M1		  ]	nuL(IBND(j))				indep
*						[M1+1		, M1+M2		  ]	nuL(IBND(j))				dep
*		j = IACT(i) =	[M1+M2+1	, 2*M1+M2	  ]	nuU(IBND(j-M1-M2))			indep
*						[2*M1+M2+1	, 2*M1+2*M2	  ]	nuU(IBND(j-M1-M2))			dep
*						[2*M1+2*M2+1, 2*M1+2*M2+M3]	lambda(dep)(j-2*M1-2*M2)
*						[	  2*M1+2*M2+M3+1	  ]	+-kapa	
*
*	UR(N+1)		O		R*8		The non-zero Lagrange multpliers for the QP constraints.
*								UR(i) is the Lagrange multplier for the active constraint
*								IACT(i).  These multipliers are categorized according to
*								IACT(i) (see IACT).
*	EXTRA		O		R*8		The slack parameter eta [0,1].  EXTRA = 0 for a feasible QP.
*	ITER		O		I		Gives the number of QP iterations used.
*	ISTATE(*)	S		I		Array for internal state data.  The length this array
*								must be at least as large as the value returned from
*								the function QPKWIKNEW_LISTATE(...)
*	LRW			W		I		Length of real workspace.  This length is returned
*								by the function QPKWIKNEW_LRW(...).
*	RW(LRW)		W		R*8		Real workspace
*
************************************************************************************************
*

      IMPLICIT NONE

*     **********************************************
*     *** Formal parameters

      INTEGER    N, M1, M2, M3, LDUINV_AUG, IBND(*), LDA, IYPY, WARM
     &           , MAX_ITER , NACTSTORE, IACTSTORE(N+1), INF, NACT
     &           , IACT(N+1), ITER, NUM_ADDS, NUM_DROPS, ISTATE(*), LRW

      DOUBLE PRECISION  GRAD(N), UINV_AUG(LDUINV_AUG,N+1), BL(*), BU(*)
     &                  , A(LDA,N), YPY(*), X(N), UR(N+1)
     &                  , EXTRA, RW(LRW), NUMPARAM(3)

*     ******************************************
*     *** External funciton types

      INTEGER QPKWIKNEW_LRW

*     ******************************************
*     *** Local variables (RAB 1/2/98, 1/13/98, 9/17/98)

      INTEGER I, J, N1, M1D, M2D, M3D
      DOUBLE PRECISION SMALL, VSMALL, VLARGE
*     *** Workspace partition indice variables
      INTEGER L_ISPARSE, L_ISTART, L_IPOINT
     &         , L_AINV, L_T1, L_T2, L_R, L_XX

*     *****************************************
*     *** Executable statements

      M1D = MAX(1, M1 + M2)
      M2D = MAX(1, M2 + M3)
      M3D = MAX(1, M1 + M2 + M3)

      SMALL  = NUMPARAM(1)
      VSMALL = NUMPARAM(2)
      VLARGE = NUMPARAM(3)

      N1 = N + 1

*     *** Reset the infeasibility

      EXTRA = 0.0

*     *** Set the iteration counter.
      ITER       = 0
      NUM_ADDS   = 0
      NUM_DROPS  = 0

*     *** Assert the length of the workspaces

      IF ( LRW .lt. QPKWIKNEW_LRW(N,M1,M2,M3) ) THEN
       INF = -2
       GOTO 10000
      ENDIF 

*     *** If the QP is unconstrained, just find the solution and return.
*     *** Replaced with BLAS, RAB (9/17/98)

      If ( M1+M2+M3.eq.0 ) Then
        NACT = 0
*       *** X = - inv(H) * GRAD = inv( U' * U ) * GRAD
*       *** X = - UINV * UINV' * GRAD
        CALL DCOPY( N, GRAD, 1, X, 1 )
        CALL DSCAL( N, -1.0, X, 1 )
        CALL DTRMV('Upper','Trans','Nonunit', N, UINV_AUG(2,1)
     &          , LDUINV_AUG, X, 1 )
        CALL DTRMV('Upper','No trans','Nonunit', N, UINV_AUG(2,1)
     &          , LDUINV_AUG, X, 1 )
        Goto 10000
      EndIf

*     *** Define the storage space

*     *** Partition ISTATE
      L_ISPARSE    = 1
      L_ISTART     = L_ISPARSE       + 1
      L_IPOINT     = L_ISTART        + M2D+1

*     *** Partition RW
      L_AINV       = 1
      L_T1         = L_AINV     + M3D+1
      L_T2         = L_T1       + N1
      L_R          = L_T2       + N1
      L_XX         = L_R        + (3*(N+1)+(N+1)*(N+1))/2

*     *** Solve the QP

      CALL QPNEW( N,M1,M2,M3,M1D,M2D,M3D,GRAD,IBND,BL,BU,A,LDA,X,NACT
     &      , IACT,UR,YPY,INF,IACTSTORE,ISTATE(L_ISTART)
     &      , ISTATE(L_IPOINT),UINV_AUG,LDUINV_AUG,RW(L_AINV),RW(L_T1)
     &      , RW(L_T2),RW(L_R),RW(L_XX),IYPY,EXTRA,MAX_ITER,ITER
     &      , NUM_ADDS,NUM_DROPS,WARM,ISTATE(L_ISPARSE),NACTSTORE
     &      , NUMPARAM )

10000 Continue

      Return
      End

*========================================================================

      SUBROUTINE QPKWIKNEW_INIT_UINV_AUG (N,BIGM,UINV,LDUINV,UINV_AUG
     &                                ,LDUINV_AUG)

*     This subroutine is just agumenting the inverse cholesky matrix of H
*     for the penalty parameter extra.  So H_AUG = [M, 0 ; 0, H] where
*     M = VLARGE.  Here:
*
*        H_AUG      = [M,         0; 0, L*L' ]
*
*        inv(H_AUG) = [1/M,       0; 0, inv(L') * inv(L)]
*
*                   = [1/sqrt(M), 0; 0, inv(L') ] * [1/sqrt(M), 0; 0, inv(L) ]
*
*                   = [1/sqrt(M), 0; 0, UINV    ] * [1/sqrt(M), 0; 0, UINV'  ]
*
*                   = UINV_AUG * UINV_AUG'
*

      IMPLICIT NONE

*     ********************************
*     *** Formal Parameters

      INTEGER N, LDUINV, LDUINV_AUG      
      DOUBLE PRECISION  BIGM, UINV(LDUINV,N), UINV_AUG(LDUINV_AUG,N+1)

*     ******************************************
*     *** Local variables

      INTEGER N1, I, J
      DOUBLE PRECISION SMALL, VSMALL, VLARGE

*     ******************************************
*     *** Executable statements

      N1 = N + 1

      UINV_AUG(1,1) = 1.0/DSQRT(BIGM)
      Do J = 2, N1
        UINV_AUG(1,J) = 0.0
      EndDo
      Do I = 2, N1
        Do J = 1, I-1
          UINV_AUG(I,J) = 0.0
        EndDo
        Do J = I, N1
          UINV_AUG(I,J) = UINV(I-1,J-1)
        EndDo
      EndDo

      RETURN
      END

*========================================================================

      INTEGER FUNCTION QPKWIKNEW_LISTATE(N,M1,M2,M3)

*     This function calculates the length of the internal state integer
*     array for QPKWIKNEW.  See the code below to determine what it is apriori

      IMPLICIT NONE

*     *** Formal parameters
      INTEGER N,M1,M2,M3

*     *** Local
      INTEGER M2D

      M2D = MAX0(1,M2+M3)

      IF ( M1 + M2 + M3 .GT. 0 ) THEN
        QPKWIKNEW_LISTATE =
*                     *** ISPARSE
     &                  1
*                     *** ISTART
     &                + M2D + 1
*                     *** IPOINT
     &                + M2D * N
      ELSE
        QPKWIKNEW_LISTATE = 1
      ENDIF

      RETURN
      END

*========================================================================

      INTEGER FUNCTION QPKWIKNEW_LRW(N,M1,M2,M3)

*     This function calculates the length of the real (double precision)
*     workspace RW for QPKWIKNEW.  See the code below to determine what it is apriori

      IMPLICIT NONE

*     *** Formal parameters
      INTEGER N,M1,M2,M3

*     *** Local
      INTEGER N1, M3D
      
      N1 = N + 1
      M3D = MAX0(1,M1+M2+M3)

      IF ( M1 + M2 + M3 .GT. 0 ) THEN
        QPKWIKNEW_LRW =
*                     *** AINV
     &                M3D + 1
*                     *** T1
     &                + N1
*                     *** T2
     &                + N1
*                     *** R
     &                + (3*(N+1)+(N+1)*(N+1))/2
*                     *** XX
     &                + N
      ELSE
        QPKWIKNEW_LRW = N
      ENDIF

      RETURN
      END

*========================================================================

*************************************
*       QPKWIKNEW Implementation    *
*************************************

*=====================================================================

      INTEGER FUNCTION JSNEW(J)

*     This function is used to access R.  Figure this out.

      IMPLICIT NONE
      INTEGER J
      JSNEW = ((J-1)+(J-1)*(J-1))/2
      RETURN
      END

*========================================================================

      SUBROUTINE QPNEW (N,M1,M2,M3,M1D,M2D,M3D,GRAD,IBND,BL,BU,A,LDA,X
     &           , NACT,IACT,UR,YPY,INF,IACTSTORE,ISTART,IPOINT
     &           , Z,LDZ,AINV,T1,T2,R,XX,IYPY,EXTRA,MAX_ITER,ITER
     &           , NUM_ADDS,NUM_DROPS,WARM,ISPARSE,NACTSTORE,NUMPARAM )

*     *** This subroutine implements the QPKWIKNEW algorithm.

      IMPLICIT NONE

*     ******************************************
*     *** Formal parameters

      INTEGER N, M1, M2, M3, M1D, M2D, M3D, IBND(M1D), LDA, NACT
     &  , IACT(N+1), INF, IACTSTORE(N+1), ISTART(M2D+1)
     &  , IPOINT(M2D*N), LDZ, IYPY, MAX_ITER, ITER, NUM_ADDS
     &  , NUM_DROPS, WARM, ISPARSE, NACTSTORE

      DOUBLE PRECISION GRAD(N), BL(M1D), BU(M3D), A(LDA,N), X(N)
     &  , UR(N+1), YPY(M1D), Z(LDZ,N+1), AINV(M3D+1), T1(N+1)
     &  , T2(N+1), R((3*(N+1)+(N+1)*(N+1))/2), XX(N), EXTRA
     &  , NUMPARAM(3)

*     ******************************************
*     *** External procedure types

      INTEGER JSNEW

*     ***************************************************************
*     *** Local variables(ordered as they first appear)

      INTEGER           N1
      DOUBLE PRECISION  SUMY, SMALL, VSMALL, VLARGE
      INTEGER           M12, M23, M123
      INTEGER           ICHECK
      INTEGER           I, J, II
      DOUBLE PRECISION  SUM
      INTEGER           KDROP, IFLAG, KSTART
      DOUBLE PRECISION  SUMNORM, CVMAX, RES
      INTEGER           KNEXT, IFINISH, IBEGIN
      DOUBLE PRECISION  TEMP
      INTEGER           INDEX
      DOUBLE PRECISION  PARNEW
      INTEGER           LFLAG
      DOUBLE PRECISION  SUMA, SUMB, SUMC, TEMPA, TEMPB
      INTEGER           IKNEXT, JJ, JN
      DOUBLE PRECISION  PARINC, STEP, RATIO
      INTEGER           ICOUNT
      DOUBLE PRECISION  XMIN, BOTTOM
      INTEGER           IWARM, ITEMP, ITEMPP

*     ***************************************************************
*     *** Executable statements

*     ************************************
*     *** [0] Initialize the algorithm

*     *** Added by RAB 1/13/98
      SMALL  = NUMPARAM(1)
      VSMALL = NUMPARAM(2)
      VLARGE = NUMPARAM(3)

*     *** Added by RAB 1/5/98 to replace common block variable
      N1 = N + 1

*     *** Check for sparsity of constraints

      M12 = M1 + M2
      M23 = M2 + M3
      M123 = M12 + M3

      CALL QPKWIK_PRINT_INPUT( N, M1, M2, M3, M1D, M2D, M3D, GRAD
     &					, Z, LDZ, IBND, BL, BU, A, LDA, YPY
     &					, INF, SMALL, VSMALL, VLARGE, N1
     &					, M12, M23, M123 )

*     *****************************************************************
*     *** [0.1] Initialize sparse access to A
      If ( M23.ne.0 ) THEN
          CALL INITNEW (N,M2D,A,ISTART,IPOINT,INF,ISPARSE)
          CALL QPKWIK_PRINT_SPARSITY(N,M2D,ISPARSE,ISTART,IPOINT)
      EndIf

*     *** ICHECK is a flag for if we should check for violation of
*     *** variable bounds and inequalities when looking for the most
*     *** violated constriant to add the the active set.  If ICHECK.eq.1
*     *** then we check all the constraints, if ICHECK.eq.0 then we only check
*     *** the equality constraints.  If there are equality constraints then
*     *** start out by only checking them first.  Once all of the equality
*     *** constraints have been added then ICHECK will be set to 1 and
*     *** then the variable bounds and inequality constriants will be
*     *** checked.
      ICHECK = 0
      If ( M12.eq.M123 ) ICHECK = 1

*     *****************************************************************
*     *** [0.2] Find the reciprocals of the constraint normal, including the simple
*     *** bounds.  Stop if a constraint is infeasible due to a zero normal.
*     *** The constraint normals are used to scale the constraint values
*     *** in determining which constraint is most violated.  In the parial
*     *** QP solver, these do not need to be calculated.

*     *** When AINV(I) is set to 0.0 we are saying that we will concider that
*     *** the constraint is always satisfied.  This happens whenever we have
*     *** an equality constraint with a near zero norm and the rhs is zero.

      Do I = 1, M1
        AINV(I) = 1.0
      EndDo
      Do I = 1, M23
        II = M1 + I
        SUM = 0.0
        Do J = ISTART(I), ISTART(I+1)-1
*         *** Need to calculate the norm of row I of A
          SUM = SUM + A(I,IPOINT(J))**2
        EndDo
        SUM = DSQRT(SUM)
        If ( SUM.le.SMALL ) Then
          If ( I.le.M2 ) Then
            If ( BL(II).le.SMALL .and.BU(II).ge.-SMALL ) Then
*             *** This is an inequality constriant that is esentially a
*             *** equality constraint with the rhs =0, BL = BU = 0.
              AINV(II) = 0.0
              GOTO 10
            EndIf
          Else
            If ( DABS(BU(II)) .le. SMALL ) Then
*             *** This is an equality constriant with a near zero norm
*             *** and a near zero rhs so we will ignore this constraint.
              AINV(II) = 0.0
              GOTO 10
            EndIf
          EndIf
        EndIf
        If ( SUM.eq.0 ) Then
          AINV(II) = 0
          GOTO 10
        EndIf
        AINV(II) = 1.0/SUM
 10     CONTINUE
      EndDo
*     *** RELAX_SPECIAL
*     *** Constraint for EXTRA >= 0
      AINV(M123+1) = 1.0

*     *****************************************************************
*     *** [0.3] Initialize to determine the unconstrained optimum.

*     *** RELAX_SPECIAL

*     *** Set EXTRA equal to zero by making the constraint EXTRA>=0 active.
*     *** This section will also be used for iterative refinement i.e.
*     *** if X changes too drastically.

      NACT = 1
      IACT(1) = M12 + M123 + 1
      UR(1) = VLARGE
      AINV(M123+1) = -1.0
      R(1) = 1.0/DSQRT(VLARGE)

*     *****************************************************************
*     *** (20)  I am not quite sure what 20 does but it updates X and UR
*     *** some how.


 20   Do I = 1, N
        X(I) = 0.0
      EndDo

*     *** RELAX_SPECIAL

      EXTRA = 0.0

*     *** RELAX_SPECIAL

*     *** If we are finding the unconstrained minimum for the
*     *** first iteration then goto it.  Here we are using
*     *** as a flag the fact that the relaxation constraint
*     *** is the only active constriant.
      If ( NACT.eq.1 .and. IACT(1).eq.M12+M123+1 ) GOTO 25

*     *** I am not sure what this code is doing.

      Do I = 1, NACT
        II = IACT(I)
        If ( II.le.M12 ) Then
*         *** For lower variable bounds and inequalities
          T1(I) = BL(II)
        ElseIf ( II.le.M12+M123 ) Then
*         *** For upper variable bounds and inequalities and equalities
          T1(I) = -BU(II-M12)
        Else
*         *** RELAX_SPECIAL
*         *** For the constriant EXTRA >= 0.0
          T1(I) = 0.0
        EndIf
      EndDo

      Do I = 1, NACT
        II = JSNEW(I)
        SUM = 0.0
        Do J = 1, I-1
          SUM = SUM + R(II+J)*T2(J)
        EndDo
        T2(I) = (T1(I)-SUM) / R(II+I)
      EndDo

*     *** Update X
      Do I = 1, NACT
*       *** RELAX_SPECIAL
        EXTRA = EXTRA + Z(1,I)*T2(I)
        Do J = 1, N
          X(J) = X(J) + Z(1+J,I)*T2(I)
        EndDo
      EndDo

*     *** Update UR
      Do I = 1, NACT
        SUM = 0.0
        Do J = 1, N
          SUM = SUM + Z(J+1,I)*GRAD(J)
        EndDo
*       *** RELAX_SPECIAL
        UR(I) = SUM + VLARGE*Z(1,I) + T2(I)
      EndDo

      Do I = NACT, 1, -1
        SUM = 0.0
        Do J = I+1, NACT
          SUM = SUM + R(I+JSNEW(J))*UR(J)
        EndDo
        UR(I) = (UR(I)-SUM) / R(I+JSNEW(I))
      EndDo

*     *** (25) RAB 3/6/98: I don't know what this does for 
*     *** an intermediate iteration.
*     ***
*     *** For the first iteration (25) will calculate the
*     *** unconstrained minimum
*     ***     X_AUG = - inv(H_AUG) * GRAD_AUG
*     ***           = - Z * Z' * GRAD_AUG
*     ***
*     ***     where X_AUG = [X ; EXTRA]
*     ***           GRAD_AUG = [GRAD ; VLARGE]

 25   CONTINUE

*     *** Compute T2 = Z * GRAD_AUG'
*     *** where GRAD_AUG = [GRAD ; VLARGE]
      Do I = NACT+1, N1
        SUM = 0.0
        Do J = 1, N
          SUM = SUM + Z(1+J,I)*GRAD(J)
        EndDo
*       *** RELAX_SPECIAL
        T2(I) = SUM + Z(1,I)*VLARGE
      EndDo
      SUM = 0.0

*     *** Compute X_AUG = Z' * T2

*     *** RELAX_SPECIAL
      Do J = NACT+1, N1
        SUM = SUM + Z(1,J)*T2(J)
      EndDo
      EXTRA = EXTRA - SUM

      Do I = 1, N
        SUM = 0.0
        Do J = NACT+1, N1
          SUM = SUM + Z(1+I,J)*T2(J)
        EndDo
        X(I) = X(I) - SUM
      EndDo

*     *** For each active variable bound, set it to its active bound.
      Do I = 1, NACT
        II = IACT(I)
        If ( II.le.M1 ) X(IBND(II)) = BL(II)
        II = II - M12
        If ( II.gt.0 .and. II.le.M1 ) X(IBND(II)) = BU(II)
      EndDo

      Do I = 1, N
        XX(I) = X(I)
      EndDo

*     *** If an inequality constraint has a negative multiplier, delete it.
      Do I = 1, NACT
        II = IACT(I)
        If ( II.le.2*M12 .and. UR(I).lt.0.0 ) Then
          KDROP = I
          IFLAG = 1
          If ( ITER .eq. MAX_ITER ) Then
             INF = -3
             Goto 1000
          EndIf
		  ITER = ITER + 1
          NUM_DROPS = NUM_DROPS + 1
          CALL DROPNEW (N,M12,M3D,NACT,IACT,AINV,KDROP,R,Z,UR
     &      ,IFLAG, NUMPARAM)
          GOTO 20
        EndIf
      EndDo

*      *** RAB: 5/11/01, We don't need this now
*      If ( INF.eq.0 ) NACTSTORE = 0

*     *** Added by RAB 1/5/98 to replace the warm start condition ALFA == 1.0
      If ( NACTSTORE.gt.0 .and. WARM.ne.1 ) NACTSTORE = 0

*     *** KSTART is a counter for which constraint we are next adding to
*     *** the active set.
      KSTART = 0

*     *** (30) Pick a constraint to add to the active set.
*     *** For the parital QP solver I need to add a callback
*     *** function to allow an external to entity determine which
*     *** constraint to add.  For now:
*     *** Look for the greatest constraint violation or
*     *** get the next constraint from the previous active set
*     *** for the warm start.

 30   CONTINUE

      CALL QPKWIK_PRINT_ITERATION_INFO( 30, N, M1, M2, M3, M1D, M2D
     &	, M3D, X, NACT, IACT, UR, IACTSTORE, Z, LDZ
     &  , AINV, T1, T2, R, XX, IYPY, EXTRA, WARM, NACTSTORE, SUMY
     &  , ICHECK, I, J, II, SUM, KDROP, IFLAG, KSTART, SUMNORM, CVMAX
     &  , RES, KNEXT, IFINISH, IBEGIN, TEMP, INDEX, PARNEW, LFLAG, SUMA
     &  , SUMB, SUMC, TEMPA, TEMPB, IKNEXT, JJ, JN, PARINC, STEP
     &  , RATIO, ICOUNT, XMIN, BOTTOM, IWARM, ITEMP, ITEMPP )

*     *** CVMAX records the greatest scaled constraint violation.
      CVMAX = VSMALL

*     *** RELAX_SPECIAL
*     *** If we made the constraint EXTRA >= 0.0 inactive to allow
*     *** for a relaxation of the feasible region then this will
*     *** always be the most violated constraint so lets pick it.
      If ( AINV(M123+1).gt.0.0 .and. -EXTRA.gt.CVMAX ) Then
        CVMAX = -EXTRA
        RES = EXTRA
        KNEXT = M12 + M123 + 1
        GOTO 40
      EndIf

*     *** RAB 3/5/98 : It looks like this part mixes the warm start
*     *** option with the search for the most violated constraint in
*     *** order to pick the next constraint to add to the active set.
*     *** Could they have been more confusing?  Is this a case of code
*     *** optimization that has made the code unmaintainable?

*     *** RAB 3/6/98 : The following code implements the warm start so it is important that
*     *** I understand it.

      If ( NACTSTORE.eq.0 ) Then
*       *** If ICHECK.eq.0 then only check the equality constraints
        If ( ICHECK.eq.0 ) IBEGIN = M12+1
*       *** If ICHECK.eq.1 then check variable bounds, equality and inequality constraints.
        If ( ICHECK.eq.1 ) IBEGIN = 1
        IFINISH = M123
      Else
*       *** Pick the next memorized constriant to check.
        IBEGIN = KSTART+1
        IFINISH = NACTSTORE
      EndIf

      Do I = IBEGIN, IFINISH
        If ( NACTSTORE.eq.0 ) Then
          II = I
        Else
          II = IACTSTORE(I)
*         *** We end up checking both upper and lower bounds so we just
*         *** need to know which constraint to check (upper and lower bounds).
          If ( II.gt.M12 ) II = II - M12
        EndIf
*       *** RELAX_SPECIAL
*       *** AINV(II) <= 0 is a flag to skip that constriant. (such as for EXTRA >= 0?). 
		If ( AINV(II).le.0.0 ) GOTO 35

*       *** Compute c(x) = x or A(i,:) * x.
        If ( II.le.M1 ) Then
          TEMP = X(IBND(II))
        Else
          TEMP = 0.0
*         *** Calculate the product of row II-M1 of A by the vector X
          Do J = ISTART(II-M1), ISTART(II-M1+1)-1
            TEMP = TEMP + A(II-M1,IPOINT(J))*X(IPOINT(J))
          EndDo
        EndIf

*       *** Compute constriant violations.

        If ( II.le.M12 ) Then
*         *** RELAX_SPECIAL
*         *** Check the lower bound on a variable or inequality.
          If ( IYPY.eq.0 ) YPY(II) = -DABS(YPY(II))
          SUM = TEMP - BL(II) - YPY(II)*EXTRA
          INDEX = II
        EndIf
        If ( II.gt.M12 .or. SUM.gt.0.0 ) Then
          If ( II.le.M12 ) Then
*           *** RELAX_SPECIAL
*           *** The lower bound was not active so check the upper bound
*           *** on the inequality.
            If ( IYPY.eq.0 ) YPY(II) = DABS(YPY(II))
            SUM = BU(II) - TEMP + YPY(II)*EXTRA
          Else
*           *** This is an equality constraint and we do not concider the
*           *** the relaxation?  This is technically not correct, is it?
            SUM = BU(II) - TEMP
          EndIf
          INDEX = M12 + II
          If ( II.le.M12 .and. SUM.gt.0.0 ) GOTO 35
        EndIf
*       *** Make negative inequality violation positive so we can
*       *** compare with an equality violation.  Also scale it
*       *** with the norm of the constraint gradient.
        SUMNORM = -SUM*AINV(II)
*       *** If this is an equality constraint, then any violation needs
*       *** to be concidered.  In this case, equality constraints are
*       *** indeed handled as double bounded inequalities.
        If ( II.gt.M12 ) SUMNORM = DABS(SUMNORM)

        If ( SUMNORM.le.CVMAX ) GOTO 35
*       *** Mark this as the next constriant to add to the active set.
        CVMAX = SUMNORM
        RES = SUM
        KNEXT = INDEX
        If ( NACTSTORE.gt.0 ) Then
*         *** If we are doing a warm start then increment to
*         *** the next memorized constraint to check the next time
*         *** a constraint is to be added to the active set.
          KSTART = I
          GOTO 40
        EndIf
        If ( II.gt.M12 ) GOTO 40
 35     CONTINUE
      EndDo

      If ( NACTSTORE.gt.0 ) Then
*       *** If you get here then you have looked through all of the
*       *** memorized constriants so just switch to the normal mode
*       *** and select the most violated constraint.
        NACTSTORE = 0
        GOTO 30
      EndIf

      If ( ICHECK.eq.0 .and. M12.ne.0 ) Then
*       *** If you get here then all of the equality constraints
*       *** are satisfied so start checking the variable bounds
*       *** and inequality constraints.
        ICHECK = 1
        GOTO 30
      EndIf

*     *** Test for convergence.  If the maximum violation is very small, then
*     *** this is the correct active set and we are finished.

      If ( CVMAX.le.VSMALL ) GOTO 700

*     *** (40) Store the constraint normal for KNEXT in T1

 40   CONTINUE

      CALL QPKWIK_PRINT_ITERATION_INFO( 40, N, M1, M2, M3, M1D, M2D
     &	, M3D, X, NACT, IACT, UR, IACTSTORE, Z, LDZ
     &  , AINV, T1, T2, R, XX, IYPY, EXTRA, WARM, NACTSTORE, SUMY
     &  , ICHECK, I, J, II, SUM, KDROP, IFLAG, KSTART, SUMNORM, CVMAX
     &  , RES, KNEXT, IFINISH, IBEGIN, TEMP, INDEX, PARNEW, LFLAG, SUMA
     &  , SUMB, SUMC, TEMPA, TEMPB, IKNEXT, JJ, JN, PARINC, STEP
     &  , RATIO, ICOUNT, XMIN, BOTTOM, IWARM, ITEMP, ITEMPP )

      Do I = 1, N1
        T1(I) = 0.0
      EndDo
      If ( KNEXT.eq.M12+M123+1 ) Then
*       *** RELAX_SPECIAL (T(1) is for the EXTRA >= 0 constraint)
        T1(1) = 1.0
      Else
        If ( KNEXT.le.M12 ) IKNEXT = KNEXT
        If ( KNEXT.gt.M12 ) IKNEXT = KNEXT - M12
        If ( IKNEXT.le.M1 ) Then
          T1(1+IBND(IKNEXT)) = 1.0
        Else
          Do I = ISTART(IKNEXT-M1), ISTART(IKNEXT-M1+1)-1
*           *** Need to assign a vector to row IKNEXT-M1 of A
            T1(1+IPOINT(I)) = A(IKNEXT-M1,IPOINT(I))
          EndDo
        EndIf
*       *** RELAX_SPECIAL
        If ( IKNEXT.le.M12 ) T1(1) = -YPY(IKNEXT)
*       *** This shows that the standard from for the inequality constriants
*       *** are A*x - bl >=0 and bu - A*x >= 0 ? 
        If ( KNEXT.gt.M12 ) Then
          Do I = 1, N1
            T1(I) = -T1(I)
          EndDo
        EndIf
      EndIf

*     *** Form the scalar product of the constraint normal for KNEXT with each column
*     *** of Z and store this as the last column of R.  PARNEW will become the
*     *** Lagrange multiplier of the new constraint.
*     *** RAB 3/9/98: This should be the BLAS call DGEMV(...).

      JN = JSNEW(NACT+1)
      Do I = 1, N1
        R(JN+I) = 0.0
        Do J = 1, N1
          R(JN+I) = R(JN+I) + Z(J,I)*T1(J)
        EndDo
      EndDo

      PARNEW = 0.0

*     *** Apply Givens rotations to make positions NACT+2 to N of R equal to zero
*     *** so that R is a (NACT+1)*(NACT+1) upper triangular matrix.

      If ( NACT.eq.N1 ) Then
        LFLAG = 4
        GOTO 50
      EndIf
      If ( ITER .eq. MAX_ITER ) Then
         INF = -3
         Goto 1000
      EndIf
	  ITER = ITER + 1
      NUM_ADDS = NUM_ADDS + 1
      CALL GIVENSNEW (N,NACT,R,Z) 

*     *** Figure out if its necessary to delete a constraint.

*     *** RAB 3/6/98 : What are the meanings of LFLAG.
*     *** I don't know what is going on here?

      If ( NACT.eq.0 .or. NACTSTORE.gt.0 ) GOTO 70
      
      SUMA = 0.0
      SUMB = 0.0
      SUMC = 0.0
      
      Do I = 1, N1
        SUMA = SUMA + T1(I)*Z(I,NACT+1)
        SUMB = SUMB + DABS(T1(I)*Z(I,NACT+1))
        SUMC = SUMC + DABS(Z(I,NACT+1))
      EndDo
      TEMPA = SUMB + 0.1D0*DABS(SUMA)
      TEMPB = SUMB + 0.2D0*DABS(SUMA)
      If ( TEMPA.le.SUMB .or. TEMPB.le.TEMPA ) Then
        LFLAG = 4
        GOTO 50
      EndIf
      IKNEXT = KNEXT
      If ( KNEXT.gt.M12 ) IKNEXT = IKNEXT - M12
      SUMC = SUMC / AINV(IKNEXT)
      TEMPA = SUMC + 0.1D0*DABS(SUMA)
      TEMPB = SUMC + 0.2D0*DABS(SUMA)
      If ( SUMC.lt.TEMPA .and. TEMPA.lt.TEMPB ) GOTO 70
      LFLAG = 1
      CALL LAGRANGENEW (LFLAG,N,M12,M123,NACT,T2,R,UR,RES,RATIO
     &               ,IACT,KDROP)
      Do I = 1, N1
        SUMA = T1(I)
        SUMB = DABS(SUMA)
        Do J = 1, NACT
          TEMP = 0.0
          JJ = IACT(J)
          If ( JJ.gt.M12 ) JJ = JJ - M12
          If ( JJ.eq.M123+1 .and. I.eq.1 ) TEMP = 1.0
          If ( JJ.le.M123 ) Then
            If ( JJ.le.M1 ) Then
              If ( I.eq.1 ) TEMP = YPY(JJ)
              If ( I.gt.1 .and. I-1.eq.IBND(JJ) ) TEMP = 1.0
            ElseIf ( JJ.le.M12 ) Then
              If ( I.eq.1 ) TEMP = YPY(JJ)
*                                  *** Need to extract an element of A
              If ( I.gt.1 ) TEMP = A(JJ-M1,I-1)
            Else
*                                  *** Need to extract an element of A
              If ( I.gt.1 ) TEMP = A(JJ-M1,I-1)
            EndIf
            If ( IACT(J).gt.M12 ) TEMP = -TEMP
          EndIf
          TEMP = TEMP*T2(J)
          SUMA = SUMA - TEMP
          SUMB = SUMB + DABS(TEMP)
        EndDo
        TEMPA = SUMB + 0.1D0*DABS(SUMA)
        TEMPB = SUMB + 0.2D0*DABS(SUMA)
        If ( SUMB.lt.TEMPA .and. TEMPA.lt.TEMPB ) GOTO 70
      EndDo
      LFLAG = 3

 50   CALL LAGRANGENEW (LFLAG,N,M12,M123,NACT,T2,R,UR,RES,RATIO
     &               ,IACT,KDROP)

*     *** Exit if the constraints are inconsistent (INF = -1)

      If ( KDROP.eq.0 ) Then
        VSMALL = VSMALL + DABS(CVMAX)
        If ( DABS(VSMALL).gt.VLARGE ) THEN
          INF = -1
          GOTO 1000
        EndIf
        GOTO 30
      EndIf

      PARINC = RATIO
      PARNEW = RATIO

*     *** Revise the Lagrange multipliers for the active constraints

 60   Do I = 1, NACT
        UR(I) = UR(I) - PARINC*T2(I)
        If ( IACT(I).le.2*M12 ) UR(I) = DMAX1(0.0d0,UR(I))
      EndDo
      If ( KDROP.eq.0 ) GOTO 80

*     *** Delete the constraint to be dropped and shift the scalar products

      IFLAG = 2
      If ( ITER .eq. MAX_ITER ) Then
         INF = -3
         Goto 1000
      EndIf
	  ITER = ITER + 1
      NUM_DROPS = NUM_DROPS + 1
      CALL DROPNEW (N,M12,M3D,NACT,IACT,AINV,KDROP,R,Z,UR,IFLAG
     &           ,NUMPARAM)

*     *** Calculate the changes to the Lagrange multipliers as well as the
*     *** step in the primal space (to reach the violated constraints)
*     ***     step_constraint = -res/sumy^2
*     *** and the largest step in the dual space (before one of the constraints
*     *** must be dropped)
*     ***     step_dual = min(ur(i),t2(i)), t2(i)>0
*     *** Do not calculate step_constraint immediately to prevent floating overflows

 70   CONTINUE

      CALL QPKWIK_PRINT_ITERATION_INFO( 70, N, M1, M2, M3, M1D, M2D
     &	, M3D, X, NACT, IACT, UR, IACTSTORE, Z, LDZ
     &  , AINV, T1, T2, R, XX, IYPY, EXTRA, WARM, NACTSTORE, SUMY
     &  , ICHECK, I, J, II, SUM, KDROP, IFLAG, KSTART, SUMNORM, CVMAX
     &  , RES, KNEXT, IFINISH, IBEGIN, TEMP, INDEX, PARNEW, LFLAG, SUMA
     &  , SUMB, SUMC, TEMPA, TEMPB, IKNEXT, JJ, JN, PARINC, STEP
     &  , RATIO, ICOUNT, XMIN, BOTTOM, IWARM, ITEMP, ITEMPP )

      SUMY = R(JSNEW(NACT+1)+NACT+1)

      If ( NACT.eq.0 ) Then
        STEP = -RES/SUMY
        PARINC = STEP/SUMY
      Else
        LFLAG = 2
        CALL LAGRANGENEW(LFLAG,N,M12,M123,NACT,T2,R,UR,RES,RATIO,
     &                IACT,KDROP)
        If ( KDROP.eq.0 ) Then
          STEP = -RES/SUMY
          PARINC = STEP/SUMY
        ElseIf ( UR(KDROP) * SUMY**2 .ge. T2(KDROP)*(-RES) ) Then
          KDROP = 0
          STEP = -RES/SUMY
          PARINC = STEP/SUMY
        Else
          RATIO = UR(KDROP)/T2(KDROP)
          TEMP = 1.0 - RATIO*SUMY**2/(-RES)
          STEP = RATIO*SUMY
          PARINC = RATIO
          RES = TEMP*RES
        EndIf
      EndIf

*     *** Update X and UR.  Drop a constraint if a full
*     *** step isn't taken.

      Do I = 1, N
        X(I) = X(I) + STEP * Z(1+I,NACT+1)
      EndDo
      EXTRA = EXTRA + STEP * Z(1,NACT+1)
      PARNEW = PARNEW + PARINC
      
      If ( NACT.ge.1 ) GOTO 60

*     *****************************************************************
*     *** (80) Add the new constraint to the active set

 80   CONTINUE

      CALL QPKWIK_PRINT_ITERATION_INFO( 80, N, M1, M2, M3, M1D, M2D
     &	, M3D, X, NACT, IACT, UR, IACTSTORE, Z, LDZ
     &  , AINV, T1, T2, R, XX, IYPY, EXTRA, WARM, NACTSTORE, SUMY
     &  , ICHECK, I, J, II, SUM, KDROP, IFLAG, KSTART, SUMNORM, CVMAX
     &  , RES, KNEXT, IFINISH, IBEGIN, TEMP, INDEX, PARNEW, LFLAG, SUMA
     &  , SUMB, SUMC, TEMPA, TEMPB, IKNEXT, JJ, JN, PARINC, STEP
     &  , RATIO, ICOUNT, XMIN, BOTTOM, IWARM, ITEMP, ITEMPP )

      NACT = NACT + 1
      UR(NACT) = PARNEW
      IACT(NACT) = KNEXT
      AINV(IKNEXT) = -AINV(IKNEXT)
      
*     *** Satisfy any bound constraint exactly and return to check remaining
*     *** constraints.
      Do I = 1, NACT
        II = IACT(I)
        If ( II.le.M1 ) X(IBND(II)) = BL(II)
        II = IACT(I) - M12
        If ( II.gt.0 .and. II.le.M1 ) X(IBND(II)) = BU(II)
      EndDo

      XMIN = 1.0
      Do I = 1, N
        If ( XX(I).ne.0.0 ) XMIN = DMIN1(DABS(X(I)/XX(I)),XMIN)
      EndDo
      If ( XMIN.lt.1.D-4 .and. IKNEXT.gt.M1 ) Then
        GOTO 20
      Else
        GOTO 30
      EndIf

*     *** Terminate Algorithm (700)
 700  CONTINUE

      CALL QPKWIK_PRINT_ITERATION_INFO( 700, N, M1, M2, M3, M1D, M2D
     &	, M3D, X, NACT, IACT, UR, IACTSTORE, Z, LDZ
     &  , AINV, T1, T2, R, XX, IYPY, EXTRA, WARM, NACTSTORE, SUMY
     &  , ICHECK, I, J, II, SUM, KDROP, IFLAG, KSTART, SUMNORM, CVMAX
     &  , RES, KNEXT, IFINISH, IBEGIN, TEMP, INDEX, PARNEW, LFLAG, SUMA
     &  , SUMB, SUMC, TEMPA, TEMPB, IKNEXT, JJ, JN, PARINC, STEP
     &  , RATIO, ICOUNT, XMIN, BOTTOM, IWARM, ITEMP, ITEMPP )

*     *** If this is not the very first QP solution then set up for a warm
*     *** start for the next QP solution.
      IWARM = 0
      If ( INF.ne.0 ) IWARM = 1
      If ( IWARM.eq.1 ) Then
*       *** Memorize the active set of constraints and sort them in
*       *** decending multiplier values scaled by the norm of the constraints.
*       *** Use T1 as a temporary sorting array.
*       *** NACSTORE and IACTSTORE are quantities that need to be saved between
*       *** calls.  They are state variables.
        Do I = 1, NACT
          II = IACT(I)
          If ( II.le.M12 ) Then
            T1(I) = -UR(I)/AINV(II)
          Else
            T1(I) = -UR(I)/AINV(II-M12)
          EndIf
*         ***
          If ( II.eq.M12+M123+1 ) T1(I) = 0.0
        EndDo
*       *** What kind of sort is this?  It is an O(NACT**2) algorithm.
*       *** Surly we can use a O(NACT*log(NACT)) algorithm.  This person
*       *** must not have known about the quick sort or merge sort.
        ICOUNT = 0
        Do I = 1, NACT
          If ( IACT(I).eq.M12+M123+1 ) GOTO 710
          ICOUNT = ICOUNT + 1
          BOTTOM = 0.0
          Do J = 1, NACT
            If ( T1(J).gt.BOTTOM ) Then
              BOTTOM = T1(J)
              ITEMP = IACT(J)
              ITEMPP = J
            EndIf
          EndDo
          IACTSTORE(ICOUNT) = ITEMP
          T1(ITEMPP) = 0.0
 710      CONTINUE
        EndDo
        NACTSTORE = ICOUNT
      EndIf

 1000 Return
      End

*========================================================================

      SUBROUTINE INITNEW (N,M2D,A,ISTART,IPOINT,INF,ISPARSE)

*     *** This subroutine sets up the sparse matrix access to A if
*     *** A is more than 50% sparse.
*

      IMPLICIT NONE

*     *************************************************************
*     *** Formal parameters

      INTEGER  N, M2D, ISTART(M2D+1), IPOINT(M2D*N), INF, ISPARSE
      DOUBLE PRECISION A(M2D,N)

*     ************************************************************
*     *** Local variables

      INTEGER           ICOUNT, I, J
      DOUBLE PRECISION  FRACT


*     *** ISPARSE must be a saved variable
*     *** From this logic it looks like the sparsity pattern is always concidered

      If ( INF.eq.0 ) ISPARSE = 1
      If ( ISPARSE.eq.2 ) RETURN

*     *** Check the sparsity of the constraints.  If less than 50% sparse use
*     *** dense version

      ICOUNT = 1
      Do I = 1, M2D
        ISTART(I) = ICOUNT
        Do J = 1, N
          If ( A(I,J) .ne. 0.0 .or. ISPARSE .eq. 1 ) Then
            IPOINT(ICOUNT) = J
            ICOUNT = ICOUNT + 1
          EndIf
        EndDo
      EndDo
      ISTART(M2D+1) = ICOUNT

*     *** On the first call (INF .eq. 0) ISPARSE is set to 1 so this if statement
*     *** will always be false and ISPARSE will be 2 after the first call.
      If ( ISPARSE.eq.0 ) Then
        FRACT = DBLE(ICOUNT-1)/DBLE(M2D*N)
        If ( FRACt.gt.0.5 ) ISPARSE = 1
      Else
        ISPARSE = 2
      EndIf

      RETURN
      END

*==========================================================================

      SUBROUTINE DROPNEW (N,M12,M3D,NACT,IACT,AINV,KDROP,R,Z,UR,IFLAG
     &   ,NUMPARAM)

*     *** This subroutine drops the constraint KDROP from the active set and
*     *** updates the QR factorization matrices Z and R.

      IMPLICIT NONE

*     *****************************  ToDo: Check the dim of R
*     *** Formal parameters

      INTEGER  N, M12, M3D, NACT, IACT(N+1), KDROP, IFLAG
      DOUBLE PRECISION AINV(M3D+1), R((3*(N+1)+(N+1)*(N+1))/2)
     &  , Z(N+1,N+1), UR(N+1), NUMPARAM(3)

*     *******************************************
*     *** Local parameters

      INTEGER           N1
      DOUBLE PRECISION  SMALL, VSMALL, VLARGE
      INTEGER           I, K, IEND, JK
      DOUBLE PRECISION  TEM, TEMP, SUM, GA, GB
      INTEGER           JI

*     *******************************************
*     *** External function types

      INTEGER JSNEW

*     *******************************************
*     *** Executable statements

*     *** Added by RAB 1/13/98
      SMALL  = NUMPARAM(1)
      VSMALL = NUMPARAM(2)
      VLARGE = NUMPARAM(3)

*     *** Added by RAB 1/5/98
      N1 = N + 1

*     *** Drop the constraint in position KDROP from the active set

      If ( IACT(KDROP).le.M12 ) Then
        AINV(IACT(KDROP)) = -AINV(IACT(KDROP))
      Else
        AINV(IACT(KDROP)-M12) = -AINV(IACT(KDROP)-M12)
      EndIf

      If ( IFLAG.eq.1 ) IEND = NACT-1
      If ( IFLAG.eq.2 .and. NACt.eq.N1 ) IEND = NACT-1
      If ( IFLAG.eq.2 .and. NACt.lt.N1 ) IEND = NACT

      Do K = KDROP, IEND

*       *** Calculate the elements of the next Givens rotation

        JK = JSNEW(K+1)
        TEM = DMAX1( DABS(R(JK+K)), DABS(R(JK+K+1)) )
        If ( TEM.gt.1.0 ) TEM = 1.0
        SUM = TEM*DSQRT((R(JK+K)/TEM)**2 + (R(JK+K+1)/TEM)**2)
        GA = R(JK+K)/SUM
        GB = R(JK+K+1)/SUM

*       *** Exchange the columns of R

        Do I = 1, K-1
          R(JSNEW(K)+I) = R(JK+I)
        EndDo
        R(JSNEW(K)+K) = SUM

*       *** Apply the rotations to the rows of R

        Do I = K+2, NACT+1
          JI = JSNEW(I)
          TEMP = GA*R(JI+K) + GB*R(JI+K+1)
          R(JI+K+1) = GA*R(JI+K+1) - GB*R(JI+K)
          R(JI+K) = TEMP
        EndDo

*       *** Apply Givens rotations to the columns of Z

        Do I = 1, N1
          TEMP = GA*Z(I,K) + GB*Z(I,K+1)
          Z(I,K+1) = GA*Z(I,K+1) - GB*Z(I,K)
          Z(I,K) = TEMP
        EndDo

*       *** Revise IACT and UR

        IACT(K) = IACT(K+1)
        UR(K) = UR(K+1)
      EndDo

*     *** If the number of active constraints is equal to the number of variables,
*     *** then the new constraint normal is shifted without any Givens rotations

      If ( IFLAG.eq.2 .and. NACt.eq.N1 ) Then
        Do I = 1, NACT
          R(JSNEW(NACT)+I) = R(JSNEW(NACT+1)+I)
        EndDo
      EndIf

      NACT = NACT - 1

      RETURN
      END

*========================================================================

      SUBROUTINE LAGRANGENEW (LFLAG,N,M12,M123,NACT,T2,R,UR,RES,RATIO
     &  , IACT,KDROP)

*     *** I don't really know what this subroutine does.

      IMPLICIT NONE

*     ***********************************************************
*     *** Formal parameters

      INTEGER    LFLAG, N, M12, M123, NACT, IACT(N+1), KDROP
      DOUBLE PRECISION T2(N+1), R((3*(N+1)+(N+1)*(N+1))/2), UR(N+1)
     &  , RES, RATIO

*     *******************************************
*     *** External function types

      INTEGER JSNEW

*     ***********************************************************
*     *** Local variables

      INTEGER           I, J
      DOUBLE PRECISION  SUM

*     ***********************************************************
*     *** Executable statements

      If ( LFLAG.eq.3 ) GOTO 10

*     *** Premultiply the vector stored in column NACT+1 by the inverse of R
*     *** and store in T2

      Do I = NACT, 1, -1
        SUM = 0.0
        Do J = I+1, NACT
          SUM = SUM + R(JSNEW(J)+I)*T2(J)
        EndDo
        T2(I) = (R(JSNEW(NACT+1)+I) - SUM) / R(JSNEW(I)+I)
      EndDo
      
      If ( LFLAG.eq.1 ) RETURN

*     *** Calculate the next constraint to drop

 10   KDROP = 0
      Do I = 1, NACT
        If ( IACT(I).gt.2*M12 .and. IACT(I).lt.M12+M123+1 ) GOTO 20
        If ( -RES*T2(I).le.0.0 ) GOTO 20
        If ( UR(I).lt.0.0 ) Then
          GOTO 20
        EndIf
        If ( KDROP.eq.0 ) Then
          KDROP = I
        ElseIf ( UR(I)*T2(KDROP).lt.UR(KDROP)*T2(I) ) Then
          KDROP = I
        EndIf
 20     CONTINUE
      EndDo
      If ( LFLAG.ne.2 ) Then
        RATIO = UR(KDROP)/T2(KDROP)
      EndIf
      
      RETURN
      END
      
*========================================================================

      SUBROUTINE GIVENSNEW (N,NACT,R,Z)

*     *** What does this subroutine do?

      IMPLICIT NONE

*     *****************************************************************
*     *** Formal parameters

      INTEGER           N, NACT
      DOUBLE PRECISION  R((3*(N+1)+(N+1)*(N+1))/2), Z(N+1,N+1)
 
*     *****************************************************************
*     *** Local variables

      INTEGER           N1, NN, I, J, K, JN
      DOUBLE PRECISION  TEMP, SUM, GA, GB

*     *******************************************
*     *** External function types

      INTEGER JSNEW

*     *****************************************************************
*     *** Executable statements

*     *** Added by RAB 1/5/98
      N1 = N + 1

*     *** Apply Givens rotations to modify Z so as to make R upper triangular

      NN = NACT + 1
      Do I = N1, NN+1, -1
        If ( R(JSNEW(NN)+I) .ne. 0.0 ) Then
          JN = JSNEW(NN)
          J = I - 1
 10       If ( R(JN+J).eq.0.0 .and. J.ne.NN ) Then
            J = J - 1
            GOTO 10
          EndIf
          TEMP = DMAX1(DABS(R(JN+J)),DABS(R(JN+I)))
          If ( TEMP.gt.1.0 ) TEMP = 1.0
          SUM = TEMP*DSQRT((R(JN+J)/TEMP)**2+(R(JN+I)/TEMP)**2)
          GA = R(JN+J)/SUM
          GB = R(JN+I)/SUM
          R(JN+J) = SUM
          Do K = N1, 1, -1
            TEMP = GA*Z(K,J) + GB*Z(K,I)
            Z(K,I) = GA*Z(K,I) - GB*Z(K,J)
            Z(K,J) = TEMP
          EndDo
        EndIf
      EndDo
      
      RETURN
      END

*========================================================================
