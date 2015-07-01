# Contributed by Robert C. Kirby
# The purpose of this example is twofold:
# 1.) Demonstrate the FIAT/PySundance interface.
#     We assume that the FIAT module (www.fenics.org/fiat)
#     is installed somewhere in the PYTHONPATH.
# 2.) Demonstrate a mixed method for Stokes flow, in which
#     the velocity and pressure are approximated by 
#     different spaces.  In this case, we use the Taylor-Hood
#     P_{k+1} velocity with P_k pressure, both are continuous

import setpath
import PySundance
import math, Numeric, string
from PySundance import *
import FIAT.Lagrange
from amesosSolver import solverParams

class LeftPointPredicate:
    def evalOp( self, x , y ):
        return math.fabs(x) < 1.0e-10

class RightPointPredicate:
    def evalOp( self, x , y ):
        return math.fabs(x-1.0) < 1.0e-10

class BottomPointPredicate:
    def evalOp( self, x , y ):
        return math.fabs(y) < 1.0e-10

class TopPointPredicate:
    def evalOp( self, x , y ):
        return math.fabs(y-1.0) < 1.0e-10

class PinPointPredicate:
    def evalOp( self , x , y ):
        return math.fabs( x ) < 1.0e-10 \
               and math.fabs(y) < 1.0e-10


pi = 4.0*math.atan(1.0)

def runtest(n,k):
    vecType = EpetraVectorType()
    npx = 1
    npy = getNProc()

	# creates n x n mesh
    mesher = PartitionedRectangleMesher(0.0,1.0,n,npx,\
                                        0.0,1.0,n,npy)

    mesh = mesher.getMesh()

	# gets velocity and pressure spaces
    velBasis = FIATScalarAdapter(FIAT.Lagrange.Lagrange,k+1)
    preBasis = FIATScalarAdapter(FIAT.Lagrange.Lagrange,k)

	# declare test and trial functions
    ux = UnknownFunction(velBasis , "ux")
    uy = UnknownFunction(velBasis , "uy")
    p  = UnknownFunction(preBasis , "p" )
    vx = TestFunction(velBasis , "vx")
    vy = TestFunction(velBasis , "vy")
    q  = TestFunction(preBasis , "q")
    x  = CoordExpr(0)
    y  = CoordExpr(1)
    dx = Derivative(0)
    dy = Derivative(1)
    grad = List(dx,dy)

    quad = GaussianQuadrature(2*k+2)

    interior = MaximalCellFilter()
    edges = DimensionalCellFilter(1)
    points = DimensionalCellFilter(0)
    left = edges.subset(LeftPointPredicate())
    right = edges.subset(RightPointPredicate())
    bottom = edges.subset(BottomPointPredicate())
    top = edges.subset(TopPointPredicate())
    pinpoint = points.subset(PinPointPredicate())
    

	# define a canned solution.
    uxex = sin( pi * x ) * cos( pi * y )
    uyex = -cos( pi * x ) * sin( pi * y )
    pex = sin( 2 * pi * x ) * sin( 2 * pi * y )
    dpex_dx = 2 * pi * cos( 2 * pi * x ) * sin( 2 * pi * y )
    dpex_dy = 2 * pi * sin( 2 * pi * x ) * cos( 2 * pi * y )
    rhs1 = 2 * pi * pi * sin( pi * x ) * cos( pi * y ) + dpex_dx
    rhs2 = - 2 * pi * pi * cos( pi * x ) * sin( pi * y ) + dpex_dy

	# in order to make a direct method behave well, we insert
	# a very small penalty term (eps * p * q) in the equations.
	# This seems to work wonders with the efficiency of Amesos
    eps = 1.e-11

    eqn = Integral( interior , \
                    (grad * vx) * (grad * ux) \
                    + (grad * vy) * (grad * uy) \
                    - p * (dx * vx + dy * vy) \
                    + q * (dx * ux + dy * uy)
                    - vx * rhs1 - vy * rhs2 \
                    + eps * p * q, quad )

    bc = EssentialBC( pinpoint , q*p , quad ) \
         + EssentialBC( left+right+bottom+top , \
                        vx * (ux-uxex) + vy * (uy-uyex) , quad ) 

    linSolver = buildSolver(solverParams)
    prob = LinearProblem( mesh , eqn , bc , \
                          List(vx,vy,q) , List(ux,uy,p) , \
                          vecType )

    print 'solving...'
    u0 = prob.solve( linSolver )
    print '...done solve'

	# we want a pressure with zero mean value.  This accomplishes that.
    pfoo = u0[2] - u0[2].integral(interior,mesh,quad)

    veldiff = (u0[0] - uxex)**2 + (u0[1] - uyex)**2
    prediff = (pfoo - pex)**2

    divvel = dx * u0[0] + dy * u0[1]

	# compute the l2 error in velocity and pressure and the 
	# l2 and l1 norms of the divergence.
    velerr = math.sqrt( veldiff.integral(interior,mesh,quad) )
    preerr = math.sqrt( prediff.integral(interior,mesh,quad) )
    divvel_l2 = divvel ** 2
    divvel_l1 = PySundance.fabs( divvel )
    sizedivvel_l2 = math.sqrt( divvel_l2.integral(interior,mesh,quad) )
    sizedivvel_l1 = divvel_l1.integral(interior,mesh,quad)

    return velerr,preerr,sizedivvel_l2,sizedivvel_l1

   
def linefit(data):
	sx = 0.0
	sy = 0.0
	sxy = 0.0
	sxx = 0.0
	
	for (x,y) in data:
		sx += x
		sxx += x*x
		sy += y
		sxy += x*y
	
	N = len(data)
	D=N*sxx - sx * sx
	a = (sxx*sy-sx*sxy)/D
	b = (N*sxy-sx*sy)/D
	
	return (a,b)


def main():
	sizes = Numeric.array( (8,16,32,64) )
	h = 1.0 / ( 1.0 * sizes ) 
	
	orders = Numeric.array( (1,2,3,4) )

	perr = Numeric.zeros( (len(orders) , len( sizes ) ) , "d" )
	uerr = Numeric.zeros( (len(orders) , len( sizes ) ) , "d" )
	divul1 = Numeric.zeros( (len(orders) , len( sizes ) ) , "d" )
	divul2 = Numeric.zeros( (len(orders) , len( sizes ) ) , "d" )


	# we loop over orders and meshes to demonstrate both h- and -p
	# convergence.  See below for a shell script that runs octave
	# to make plots of convergence rates for a problem with a 
	# canned solution

	for i in range( len( orders ) ):
		order = orders[i]
		for j in range( len( sizes ) ):
			n = sizes[j]
			uerr[i,j],perr[i,j],divul1[i,j],divul2[i,j]=runtest(n,order)


	for dat,filename in zip( (perr,uerr,divul1,divul2) , \
							 ("perr.dat","uerr.dat","divul1.dat","divul2.dat") ):
		datafile = open( filename , "w" )
		for i in range( len( orders ) ):
			print >>datafile , string.join( map(str , dat[i]) )
		datafile.close()

	datafile = open( "orders.dat" , "w" )
	print >>datafile , string.join( map( str , orders ) , "\n" )
	datafile.close()
	datafile = open( "h.dat" , "w" )
	print >>datafile , string.join( map( str , h ) , "\n" )

# main shell script
	octavemain = open( "TH.sh" , "w" )
	print >>octavemain, """
#!/bin/sh
octave perrhplot.m
octave perrkplot.m
octave uerrhplot.m
octave uerrkplot.m
octave divul1hplot.m
octave divul1kplot.m
octave divul2hplot.m
octave divul2kplot.m"""
	octavemain.close()

	# make plots of pressure error versus h, one line for each degree
	# I need to make the titles; each of these are the order of the pressure
	legends = map(lambda x: "\";P%s;\"" % (x,) , orders )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "perr(:,%s)'" % (i+1,) for i in range(len(orders)) ] , \
								  " ; " ) , )
	octavefile = open( "perrhplot.m" , "w" )
	print >>octavefile, \
"""load "h.dat";
load "perr.dat";
loglog( h' , %s , %s );
gset grid;
gset xlabel "h";
gset ylabel "||p-ph||_L2";
gset terminal postscript;
gset output "perrh.ps";
""" % (ys,legendstr)
	octavefile.close()

	# plots of pressure error versus k, one line for each mesh
	legends = map( lambda x: "\";h=%s;\"" % (x,) , h )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "perr(%s,:)" % (i+1,) for i in range(len(h)) ] , \
								  " ; " ) , )

	octavefile=open("perrkplot.m" , "w")
	print >>octavefile, \
"""load "orders.dat";
load "perr.dat";
semilogy( orders' , %s , %s );
gset grid;
gset xlabel "k";
gset ylabel "||p-ph||_L2";
gset terminal postscript;
gset output "perrk.ps";
""" % (ys,legendstr)
	octavefile.close()


	# make plots of velocity error versus h, one line for each degree
	# I need to make the titles; each of these are the order of the pressure
	legends = map(lambda x: "\";P%s;\"" % (x,) , 1+orders )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "uerr(:,%s)'" % (i+1,) for i in range(len(orders)) ] , \
								  " ; " ) , )
	octavefile=open("uerrhplot.m" , "w")
	print >>octavefile, \
"""load "h.dat";
load "uerr.dat";
loglog( h' , %s , %s );
gset grid;
gset xlabel "h";
gset ylabel "||u-uh||_L2";
gset terminal postscript;
gset output "uerrh.ps";
""" % (ys,legendstr)
	octavefile.close()

	# plots of velocity error versus k, one line for each mesh
	legends = map( lambda x: "\";h=%s;\"" % (x,) , h )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "uerr(%s,:)" % (i+1,) for i in range(len(h)) ] , \
								  " ; " ) , )

	octavefile=open("uerrkplot.m" , "w")
	print >>octavefile, \
"""load "orders.dat";
load "uerr.dat";
semilogy( (1+orders)' , %s , %s );
gset grid;
gset xlabel "k";
gset ylabel "||u-uh||_L2";
gset terminal postscript;
gset output "uerrk.ps";
""" % (ys,legendstr)
	octavefile.close()

	# make plots of l2 norm of divergence versus h, one line for each degree
	legends = map(lambda x: "\";P%s;\"" % (x,) , orders + 1)
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "divul1(:,%s)'" % (i+1,) for i in range(len(orders)) ] , \
								  " ; " ) , )
	octavefile = open("divul1hplot.m" , "w" )
	print >>octavefile, \
"""load "h.dat";
load "divul1.dat";
loglog( h' , %s , %s );
gset grid;
gset xlabel "h";
gset ylabel "||div(uh)||_L1";
gset terminal postscript;
gset output "divul1h.ps";
""" % (ys,legendstr)

	octavefile.close()

	# plots of pressure error versus k, one line for each mesh
	legends = map( lambda x: "\";h=%s;\"" % (x,) , h )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "divul1(%s,:)" % (i+1,) for i in range(len(h)) ] , \
								  " ; " ) , )

	octavefile=open("divul1kplot.m" , "w")
	print >>octavefile, \
"""load "orders.dat";
load "divul1.dat";
semilogy( (1+orders)' , %s , %s );
gset grid;
gset xlabel "k";
gset ylabel "||div(uh)||_L1";
gset terminal postscript;
gset output "divul1k.ps";
""" % (ys,legendstr)
	octavefile.close()


	# make plots of l2 norm of divergence versus h, one line for each degree
	legends = map(lambda x: "\";P%s;\"" % (x,) , orders + 1)
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "divul2(:,%s)'" % (i+1,) for i in range(len(orders)) ] , \
								  " ; " ) , )
	octavefile = open("divul2hplot.m" , "w" )
	print >>octavefile, \
"""load "h.dat";
load "divul2.dat";
loglog( h' , %s , %s );
gset grid;
gset xlabel "h";
gset ylabel "||div(uh)||_L2";
gset terminal postscript;
gset output "divul2h.ps";
""" % (ys,legendstr)

	octavefile.close()

	# plots of pressure error versus k, one line for each mesh
	legends = map( lambda x: "\";h=%s;\"" % (x,) , h )
	legendstr = "[ %s ]" % (string.join(legends,";"),)
	ys = "[ %s ]" % (string.join( [ "divul2(%s,:)" % (i+1,) for i in range(len(h)) ] , \
								  " ; " ) , )

	octavefile=open("divul2kplot.m" , "w")
	print >>octavefile, \
"""load "orders.dat";
load "divul2.dat";
semilogy( (1+orders)' , %s , %s );
gset grid;
gset xlabel "k";
gset ylabel "||div(uh)||_L2";
gset terminal postscript;
gset output "divul2k.ps";
""" % (ys,legendstr)
	octavefile.close()


if __name__ == "__main__":
    main()
    


    
