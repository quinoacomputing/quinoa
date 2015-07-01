import setpath
from PySundance import *



amesosSolverDict = {"Linear Solver" : 
                    {"Type" : "Amesos",
                     "Kernel" : "Umfpack"}
                    }

solverParams = dict2ParameterList(amesosSolverDict);
