import setpath
from PySundance import *



blockSolverDict = {"Linear Solver" :
                   {"Type" : "Block Triangular",
                    "Sub Solver" :
                    {"Linear Solver" : 
                     {"Type" : "Aztec",
                      "Method" : "GMRES",
                      "Max Iterations" : 1000,
                      "Tolerance" : 1.0e-12,
                      "Precond" : "Domain Decomposition",
                      "Subdomain Solver" : "ILU",
                      "Graph Fill" : 1,
                      "Verbosity" : 4
                      }
                     }
                    }
                   }

solverParams = dict2ParameterList(blockSolverDict);
