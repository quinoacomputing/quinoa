#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
from PySundance import *

aztecSolverDict = {"Type" : "Aztec",
                   "Method" : "GMRES",
                   "Max Iterations" : 1000,
                   "Tolerance" : 1.0e-8,
                   "Precond" : "Domain Decomposition",
                   "Subdomain Solver" : "ILU",
                   "Graph Fill" : 1,
                   "Verbosity" : 4
                   }
               

solverDict = {"NOX Solver" :
              {"Nonlinear Solver" : "Line Search Based",
               "Line Search" : {"Method" : "More'-Thuente"},
               "StatusTest"  : {"Max Iterations" : 20, "Tolerance" : 1.0e-6},
               "Linear Solver" : aztecSolverDict}}

solverParams = dict2ParameterList(solverDict);

