
// @HEADER
// ***********************************************************************
//
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ***********************************************************************
// @HEADER

#include "Didasko_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#if defined(HAVE_DIDASKO_TEUCHOS)

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ConfigDefs.hpp"

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Creating an empty parameter list looks like:
  Teuchos::ParameterList My_List;

  // Setting parameters in this list can be easily done:
  My_List.set("Max Iters", 1550);
  My_List.set("Tolerance", 1e-10);
  My_List.set("Solver", "GMRES");

  /* The templated ``set'' method should cast the input {\it value} to the
     correct data type.  However, in the case where the compiler is not casting the input
     value to the expected data type, an explicit cast can be used with the ``set'' method:
     */
  My_List.set("Tolerance", (float)(1e-10));

  /* A hierarchy of parameter lists can be constructed using {\tt Teuchos::ParameterList}.  This
     means another parameter list is a valid {\it value} in any parameter list.  To create a sublist
     in a parameter list and obtain a reference to it:
     */
  Teuchos::ParameterList& Prec_List = My_List.sublist("Preconditioner");

  // Now this parameter list can be filled with values:
  Prec_List.set("Type", "ILU");
  Prec_List.set("Drop Tolerance", 1e-3);

  // The parameter list can be queried about the existance of a parameter, sublist, or type:
  // Has a solver been chosen?
  bool solver_defined, prec_defined, tol_double, dtol_double;
  solver_defined = My_List.isParameter("Solver");
  // Has a preconditioner been chosen?
  prec_defined = My_List.isSublist("Preconditioner");
  // Has a tolerance been chosen and is it a double-precision number?
  tol_double = My_List.INVALID_TEMPLATE_QUALIFIER isType<double>("Tolerance");
  // Has a drop tolerance been chosen and is it a double-precision number?
  dtol_double = Teuchos::isParameterType<double>(Prec_List, "Drop Tolerance");

  /* The last two methods for checking the parameter type are equivalent.
     There is some question as to whether the syntax of the first type-checking
     method is acceptable to older compilers.  Thus, the second type-checking method
     is offered as a portable alternative.
     */
  // Parameters can be retrieved from the parameter list in quite a few ways:
  // Get method that creates and sets the parameter if it doesn't exist.
  int its;
  its = My_List.get("Max Iters", 1200);
  // Get method that retrieves a parameter of a particular type.
  float tol = My_List.INVALID_TEMPLATE_QUALIFIER get<float>("Tolerance");

  /* In the above example, the first ``get'' method is a safe way of
     obtaining a parameter when its existence is indefinite but required.
     The second ``get'' method should be used when the existense of the parameter
     is definite.  This method will throw an exception if the parameter doesn't exist.
     The safest way to use the second ``get'' method
     is in a try/catch block:
     */
  try {
    tol = My_List.INVALID_TEMPLATE_QUALIFIER get<float>("Tolerance");
  }
  catch (std::exception& e) {
    tol = 1e-6;
  }

  /* The second ``get'' method uses a syntax that may not be
     acceptable to older compilers.  Optionally, there is another portable templated
     ``get'' function that can be used in the place of the second ``get'' method:
     */
  try {
    tol = Teuchos::getParameter<float>(My_List, "Tolerance");
  }
  catch (std::exception& e) {
    tol = 1e-6;
  }

  // A parameter list can be sent to the output stream:
  cout<< My_List << endl;

  /* It is important to note that mispelled parameters
     (with additional space characters, capitalizations, etc.) may be ignored.
     Therefore, it is important to be aware that a given parameter has not been used.
     Unused parameters can be printed with method:
     */
  My_List.unused( cout );

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:\n"
      "--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif
