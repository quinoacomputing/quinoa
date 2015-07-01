#!/usr/bin/env python

import setpath
import PySundance

import math
from PySundance import *

import FIAT

def main():
  """Poisson example code"""
  basis = Lagrange(2)
  basis2 = FIATScalarAdapter( FIAT.Lagrange.Lagrange , 2 )

  checkbasis( basis , basis2 )
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
