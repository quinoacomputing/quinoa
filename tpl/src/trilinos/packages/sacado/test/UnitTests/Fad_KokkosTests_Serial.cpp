// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Fad_KokkosTests.hpp"

#include "Kokkos_Core.hpp"

// Instantiate tests for Serial device
using Kokkos::Serial;
VIEW_FAD_TESTS_D( Serial )

// Add a unit test verifying something from Albany compiles
TEUCHOS_UNIT_TEST(Kokkos_View_Fad, DynRankMauroDeepCopy )
{
  Kokkos::DynRankView<Sacado::Fad::DFad<double>,Kokkos::Serial> v1(
    "v1", 3, 5);
  Kokkos::DynRankView<Sacado::Fad::DFad<double>,Kokkos::LayoutRight,Kokkos::HostSpace> v2("v2", 3 , 5);

  Kokkos::deep_copy(v1, v2);

  // We're just verifying this compiles
  success = true;
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize serial
  Kokkos::Serial::initialize();
  if (!std::is_same<Kokkos::Serial, Kokkos::HostSpace::execution_space>::value)
    Kokkos::HostSpace::execution_space::initialize();

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize serial
  if (!std::is_same<Kokkos::Serial, Kokkos::HostSpace::execution_space>::value)
    Kokkos::HostSpace::execution_space::finalize();
  Kokkos::Serial::finalize();

  return res;
}
