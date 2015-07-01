#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance

import math
from PySundance import *
from random import *


def main():

  vecType = EpetraVectorType()
  ny = 10
  nx = 10
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, 1,
                                       0.0, 1.0, ny, 1);
  mesh = mesher.getMesh();

  x = CoordExpr(0);
  y = CoordExpr(1);


  interior = MaximalCellFilter()


  # Define a force field as a discrete function. For testing purposes,
  # we use a force that is linear in (x,y) so that we can do an exact
  # comparison with its interpolated values. In a real problem, the force
  # can of course be more complicated.
  # NOTE: point interpolation is currently only supported for
  # discrete functions with Lagrange(1) basis.
  L1 = Lagrange(1)
  discSpace = DiscreteSpace(mesh, BasisList(L1, L1), vecType)
  projector1 = L2Projector(discSpace, List(x+2.0*y, y-0.2*x))
  F1 = projector1.project()


  # Create a point locator object. This same object can be shared between
  # interpolator objects and sampler objects. The virtualGridSize argument
  # defines the virtual structured grid used to accelerate lookups.
  virtualGridSize = [27, 27]
  locator = AToCPointLocator(mesh, interior, virtualGridSize)

  # Create an interpolator object that will use our locator to interpolate
  # values from F1.
  forceInterpolator = CToAInterpolator(locator, F1)

  
  # Create a vector of 2D points with interleaved coordinates, i.e.,
  # points = [x_0, y_0, x_1, y_1, ..., x_N, y_N]. The "doubleVector" type
  # is a python-wrapped std::vector<double> object.
  nPts = 1500
  points = doubleVector(2*nPts)

  # Choose points from a uniform distribution on the unit square
  for i in range(nPts):
      points[2*i] = random()
      points[2*i+1] = random()

  # Create a vector for the interpolated forces.
  forces = doubleVector(2*nPts)

  # Do the interpolation
  forceInterpolator.interpolate(points, forces)

  # Compute the error in the force interpolations. Because we used
  # a linear force and are doing linear interpolation, this should
  # be very small, attributable to roundoff.
  maxError = 0.0

  for i in range(nPts):
      x0 = [points[2*i], points[2*i+1]]
      f = [x0[0] + 2.0*x0[1], x0[1] - 0.2*x0[0]]
      df = math.fabs(f[0] - forces[2*i]) + math.fabs(f[1] - forces[2*i+1])
      if df > maxError: maxError = df


  # Now change the discrete force field.
  projector2 = L2Projector(discSpace, List(x+y, x-y))
  F1 = projector2.project()

  # Tell the interpolator that the values of F1 have changed.
  # IMPORTANT: for efficiency reasons, the interpolator makes its own local
  # copy of the discrete function values in a format optimized
  # for pointwise interpolation. Therefore, it is necessary to update
  # the force interpolator any time the discrete function changes, even
  # if it has been changed "in-place" as in a nonlinear solver. Otherwise,
  # the Interpolator's local values will be out of synch with the
  # DiscreteFunction's values. Note: in principle, I can autodetect changes,
  # but haven't yet done so.
  forceInterpolator.updateField(F1)

  # Interpolate the new force field, and check errors
  forceInterpolator.interpolate(points, forces)

  for i in range(nPts):
      x0 = [points[2*i], points[2*i+1]]
      f = [x0[0] + x0[1], x0[0] - x0[1]]
      df = math.fabs(f[0] - forces[2*i]) + math.fabs(f[1] - forces[2*i+1])
      if df > maxError: maxError = df

  
  # Next we sample the density of points. First, create a sampler
  # object using the same locator defined above. A VectorType argument
  # is also necessary, because we'll be creating a new DiscreteSpace
  # of piecewise-constant functions (i.e., Lagrange(0)).
  sampler = AToCDensitySampler(locator, vecType)

  # Compute the density, weighting each particle by 1/nPts so that
  # the total mass is 1.0. This weight factor can also be used to set
  # appropriate units, e.g., converting from particles per cubic nanometer
  # to moles per liter.
  density = sampler.sample(points, 1.0/nPts)

  # Compute the total mass. This should be 1.0.
  mass = density[0].integral(interior, mesh, GaussianQuadrature(2))
  massErr = math.fabs(mass - 1.0)
  if massErr > maxError: maxError = massErr

  # Write the density field

  writer = VTKWriter("Density2D");
  writer.addMesh(mesh)
#  writer.addField("rho0", density[0])

  # ----------------------------------------------------------------------------
  # Finally, we will sample the points with weighting appropriate to axisymmetric
  # geometry, and in multiple sample steps. We will do several outer steps, resetting
  # the count to zero between each. In each outer step there are nInnerSteps inner steps,
  # during which we add into the density field.

  # Specify axis of rotation for the coordinate system
  origin = [0.0, 0.0]
  axisOfRotation = [1.0, 0.0]
  # Create a sampler for this coordinate system
  sampler = AToCDensitySampler(locator, origin, axisOfRotation, vecType)



  # Do the run over outer and inner steps
  for outerStep in range(3):
    print 'outer step=', outerStep
    # reset the counts to zero
    density2 = sampler.resetCounts()
    
    # do the inner timesteps, sampling at each step
    nInnerStep = 20
    for innerStep in range(nInnerStep) :
      # sample points uniformly from a cylinder, projecting to the (z,r) plane
      for i in range(nPts):
        points[2*i] = random()
        points[2*i+1] = math.sqrt(random())
      sampler.addToCounts(points, 1.0/(nPts*nInnerStep), density2)

    # write the current sample
    writer.addField('rhoStep%d' % outerStep, density2[0])
    

  # Write the results
  writer.write()

  # all done!  Check the error and quit.
  tol = 1.0e-10
  passFailTest(maxError, tol)




  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
