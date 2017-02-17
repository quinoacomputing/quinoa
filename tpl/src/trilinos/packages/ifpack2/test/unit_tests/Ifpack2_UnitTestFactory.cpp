// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Ifpack2_UnitTestFactory.cpp

\brief Ifpack2 Unit test for the Factory template, and some basic tests
for preconditioners it produces.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Factory.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void check_precond_basics(Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& precond, Teuchos::FancyOStream& out, bool& success)
{
  TEST_EQUALITY( precond->getNumInitialize(), 0);
  TEST_EQUALITY( precond->getNumCompute(), 0);
  TEST_EQUALITY( precond->getNumApply(), 0);

  precond->initialize();
  precond->compute();

  TEST_EQUALITY( precond->getNumInitialize(), 1);
  TEST_EQUALITY( precond->getNumCompute(), 1);
}

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Factory, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Factory factory;

  Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_ilut = factory.create("ILUT", crsmatrix);
  TEST_EQUALITY(prec_ilut != Teuchos::null, true);

  check_precond_basics(prec_ilut, out, success);

  Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_riluk = factory.create("RILUK", crsmatrix);
  TEST_EQUALITY(prec_riluk != Teuchos::null, true);

  check_precond_basics(prec_riluk, out, success);

  Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_relax = factory.create("RELAXATION", crsmatrix);
  TEST_EQUALITY(prec_relax != Teuchos::null, true);

  check_precond_basics(prec_relax, out, success);

  Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_diag = factory.create("DIAGONAL", crsmatrix);
  TEST_EQUALITY(prec_diag != Teuchos::null, true);

  check_precond_basics(prec_diag, out, success);

  Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_cheby = factory.create("CHEBYSHEV", crsmatrix);
  TEST_EQUALITY(prec_cheby != Teuchos::null, true);

  check_precond_basics(prec_cheby, out, success);
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Factory, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

}//namespace <anonymous>

