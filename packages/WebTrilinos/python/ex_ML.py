#! /bin/env python
# @HEADER
# ************************************************************************
#
#                WebTrilinos: A Web Interface to Trilinos
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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

try:
  import setpath
  import Epetra, AztecOO, ML
except:
  from PyTrilinos import Epetra, AztecOO, ML

Comm = Epetra.PyComm()

NumGlobalElements = 10
Map = Epetra.Map(NumGlobalElements, 0, Comm)
MyGlobalElements = Map.MyGlobalElements()

Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0)

for i in MyGlobalElements:
  Matrix[i, i] = 1.0

Matrix.FillComplete()

MLList = {
  "max levels"        : 3,
  "output"            : 10,
  "smoother: type"    : "symmetric Gauss-Seidel",
  "aggregation: type" : "Uncoupled"
}

Prec = ML.MultiLevelPreconditioner(Matrix, False)
Prec.SetParameterList(MLList)
Prec.ComputePreconditioner()

LHS = Epetra.Vector(Map); LHS.PutScalar(0.0)
RHS = Epetra.Vector(Map); RHS.Random()

Solver = AztecOO.AztecOO(Matrix, LHS, RHS)

Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_output, 16)
Solver.Iterate(1550, 1e-5)

del Prec

print Solver.NumIters(), Solver.TrueResidual()
