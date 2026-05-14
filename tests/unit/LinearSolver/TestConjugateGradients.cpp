// *****************************************************************************
/*!
  \file      tests/unit/LinearSolver/TestConjugateGradients.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for LinearSolver/ConjugateGradients
  \details   Unit tests for LinearSolver/ConjugateGradients
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "ConjugateGradients.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "Vector.hpp"

#include "NoWarning/cgreceiver.decl.h"
#include "NoWarning/tutsuite.decl.h"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct ConjugateGradients_common {};

//! Test group shortcuts
using ConjugateGradients_group =
  test_group< ConjugateGradients_common, MAX_TESTS_IN_GROUP >;
using ConjugateGradients_object = ConjugateGradients_group::object;

//! Define test group
static ConjugateGradients_group
  ConjugateGradients( "LinearSolver/ConjugateGradients" );

//! Charm++ chare to receive callbacks when CG tasks are complete
class CGReceiver : public CBase_CGReceiver {
  public:
    //! Constructor
    CGReceiver( const std::string& label,
                std::size_t maxit,
                tk::real tol,
                tk::real normb_ex,
                tk::real normres_ex,
                tk::CProxy_ConjugateGradients cg )
      : m_label( label ), m_maxit(maxit), m_tol(tol), m_normb_ex(normb_ex),
        m_normres_ex(normres_ex), m_cg( cg ) {}
    //! Called after CG::setup() finished
    void initialized( CkDataMsg* msg ) {
      auto normb = static_cast<tk::real*>( msg->getData() );
      received( "init, ch" + std::to_string(thisIndex),
                "normb", *normb, m_normb_ex );
      m_cg[thisIndex].solve( m_maxit, m_tol,
        CkCallback(CkIndex_CGReceiver::solved(nullptr),thisProxy[thisIndex]) );
    }
    //! Called after CG solve() finished
    void solved( CkDataMsg* msg ) {
      auto residual = static_cast<tk::real*>( msg->getData() );
      received( "solve, ch" + std::to_string(thisIndex),
                "residual", *residual, m_normres_ex );
    }
  private:
    //! Test label this is a receiver for
    std::string m_label;
    //! Max iteration count
    std::size_t m_maxit;
    //! Stop tolerance
    tk::real m_tol;
    //! Expected norm of the right hand side (b) vector
    tk::real m_normb_ex;
    //! Expected norm of the final residual
    tk::real m_normres_ex;
    //! CG solver proxy
    tk::CProxy_ConjugateGradients m_cg;
    //! Function creating a TUT test on completing a CG task
    void received( const std::string& msg,
                   const std::string& in,
                   tk::real computed,
                   tk::real expected )
    {
      // Create test result struct
      tut::test_result tr( "LinearSolver/ConjugateGradients", 1,
                           m_label + " " + msg,
                           tut::test_result::result_type::ok );
      try {
        // Evaluate test
        ensure_equals( "CG error in " + in, computed, expected, 1.0e-12 );
      } catch ( const failure& ex ) {
        tr.result = ex.result();
        tr.exception_typeid = ex.type();
        tr.message = ex.what();
      }
      // Send back a new test result, signaling the second part
      unittest::g_suiteProxy.evaluate(
        { tr.group, tr.name, std::to_string(tr.result), tr.message,
          tr.exception_typeid } );
    }
};

//! Test definitions for group

//! Test simple CG solve with DOF=1 of Laplacian in serial
template<> template<>
void ConjugateGradients_object::test< 1 >() {
  // This test spawns a new Charm++ chare which creates two new TUT tests. Each
  // test has similar test names but differ at the end, corresponding to a
  // stage of the completion of the multiple tests as part of the series. The
  // spawed tests create new test results, sending them back to the suite if
  // successful.
  set_test_name( "Laplacian 1DOF 1PE setup" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel {
    3, 13, 8, 14,
    12, 3, 13, 8,
    8, 3, 14, 11,
    12, 3, 8, 11,
    1, 2, 3, 13,
    6, 13, 7, 8,
    5, 9, 14, 11,
    5, 1, 3, 14,
    10, 4, 12, 11,
    2, 6, 12, 13,
    8, 7, 9, 14,
    13, 1, 7, 14,
    5, 3, 4, 11,
    6, 10, 12, 8,
    3, 2, 4, 12,
    10, 8, 9, 11,
    3, 1, 13, 14,
    13, 7, 8, 14,
    6, 12, 13, 8,
    9, 8, 14, 11,
    3, 5, 14, 11,
    4, 3, 12, 11,
    3, 2, 12, 13,
    10, 12, 8, 11 };

  // Mesh node coordinates for simple tet mesh above
  std::array< std::vector< tk::real >, 3 > coord {{
    {{ -0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0 }},
    {{ 0.5, 0.5, 0, -0.5, -0.5, 0.5, 0.5, 0, -0.5, -0.5, -0.5, 0, 0.5, 0 }},
    {{ -0.5, 0.5, 0, 0.5, -0.5, 0.5, -0.5, 0, -0.5, 0.5, 0, 0.5, 0, -0.5 }} }};

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );
  // Query number of nodes in mesh
  auto npoin = psup.second.size()-1;

  // Create CSR matrix based on mesh with psup
  tk::CSR A( /* DOF = */ 1, psup );

  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];

  // fill matrix with Laplacian
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );

    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    for (std::size_t a=0; a<4; ++a)
      for (std::size_t k=0; k<3; ++k)
         for (std::size_t b=0; b<4; ++b)
           A(N[a],N[b]) += J/6 * grad[a][k] * grad[b][k];
  }

  // Create RHS and solution/unknown vectors
  std::vector< tk::real > b(npoin,1.0), x(npoin,0.0);

  // Grab a node (id=0) as Dirichlet BC
  A.dirichlet( 0 );

  // Create CG solver (chare array with a single element)
  tk::CProxy_ConjugateGradients cg =
    tk::CProxy_ConjugateGradients::ckNew( A, x, b, {}, {}, {}, 1 );

  // Create receiver chare whose callbacks are called when a CG task is done
  CProxy_CGReceiver host =
    CProxy_CGReceiver::ckNew( "Laplacian 1DOF 1PE",
                              /* maxit = */ 200,
                              /* tol = * */ 1.0e-16,
                              /* normb_ex = */ 3.7416573867739413,
                              /* normres_ex = */ 3.8368230753371624e-15,
                              cg, 1 );

  // Initialize CG solve
  cg[0].setup( CkCallback(CkIndex_CGReceiver::initialized(nullptr), host[0]) );
}

//! Test simple CG solve with DOF=1 of Laplacian on 2 PEs
template<> template<>
void ConjugateGradients_object::test< 2 >() {
  // This test spawns a new Charm++ chare array of two elements which each
  // creates two new TUT tests, thus resulting four newly spawned TUT tests.
  // Each test has similar test names but differ, corresponding to a stage of
  // the completion of the multiple tests as part of the series and
  // corresponding to which chare array element they complete. The spawed tests
  // create new test results, sending them back to the suite if successful.
  set_test_name( "Laplacian 1DOF 2PE setup" );

  // Initialize psup for PE0 and PE1
  std::vector<
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > > psup {
    { { 0,1,2,4,8,9,0,2,3,7,8,0,1,3,4,5,6,7,8,9,1,2,4,6,7,0,2,3,6,9,2,6,7,8,9,2,
        3,4,5,7,9,1,2,3,5,6,8,0,1,2,5,7,9,0,2,4,5,6,8 },
    { 0,5,10,19,24,29,34,40,46,52,58 } },
    { { 0,5,11,12,4,10,11,8,9,10,7,9,12,1,5,6,8,10,11,0,4,6,7,11,12,4,5,7,8,9,10,
      11,12,3,5,6,8,9,12,2,4,6,7,9,10,2,3,6,7,8,10,12,1,2,4,6,8,9,11,0,1,4,5,6,
      10,12,0,3,5,6,7,9,11 },
    { 0,3,6,9,12,18,24,32,38,44,51,58,65,72 } } };

  // Initialize node coordinates for PE0 and PE1
  std::vector< std::array< std::vector< tk::real >, 3 > > coord {{
    {{ {{ -0.5,-0.5,-0.5,-0.5,-0.5,0.5,0,0,0,0 }},
       {{ 0.5,0.5,0,-0.5,-0.5,0,-0.5,0,0.5,0 }},
       {{ -0.5,0.5,0,0.5,-0.5,0,0,0.5,0,-0.5 }} }},
    {{ {{ -0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0 }},
       {{ 0.5,0.5,-0.5,-0.5,0.5,0.5,0,-0.5,-0.5,-0.5,0,0.5,0 }},
       {{ -0.5,0.5,0.5,-0.5,0.5,-0.5,0,-0.5,0.5,0,0.5,0,-0.5 }} }} }};

  // Initialize mesh connectivity for PE0 and PE1
  std::vector< std::vector< std::size_t > > inpoel {
    { 2,8,5,9,7,2,8,5,5,2,9,6,7,2,5,6,0,1,2,8,4,0,2,9,4,2,3,6,2,1,3,7,
      2,0,8,9,2,4,9,6,3,2,7,6,2,1,7,8 },
    { 4,8,10,6,8,6,7,9,11,5,6,12,4,10,11,6,7,6,12,9,8,10,6,9,4,11,5,6,
      3,7,12,9,8,2,10,9,1,4,10,11,6,5,7,12,11,0,5,12 } };

  // Initialize global mesh node IDs for PE0 and PE1
  std::vector< std::vector< std::size_t > > gid {
    { 0,1,2,3,4,7,10,11,12,13 },
    { 0,1,3,4,5,6,7,8,9,10,11,12,13 } };

  // Initialize local->global mesh node ID map for PE0 and PE1
  std::vector< std::unordered_map< std::size_t, std::size_t > > lid{
    { {0,0},{11,7},{1,1},{12,8},{2,2},{13,9},{3,3},{4,4},{7,5},{10,6} },
    { {0,0},{11,10},{5,4},{1,1},{3,2},{4,3},{6,5},{7,6},{8,7},{9,8},
      {10,9},{12,11},{13,12} } };

  // Initialize node communication map for PE0 and PE1
  std::vector< tk::NodeCommMap > nodecommap;
  nodecommap.push_back({});
  nodecommap.back()[1].insert( { 7,4,10,12,1,13,3,0,11 } );
  nodecommap.push_back({});
  nodecommap.back()[0].insert( { 7,4,3,10,1,12,11,0,13 } );

  // Create CG solver (empty chare array, will use dynamic insertion)
  tk::CProxy_ConjugateGradients cg = tk::CProxy_ConjugateGradients::ckNew();

  // Create receiver chare array (2 elements) whose callbacks are called when a
  // CG task is done on a PE
  CProxy_CGReceiver host =
    CProxy_CGReceiver::ckNew( "Laplacian 1DOF 2PE",
                              /* maxit = */ 100,
                              /* tol = * */ 1.0e-16,
                              /* normb_ex = */ 3.7416573867739413,
                              /* normres_ex = */ 3.8368230753371624e-15,
                              cg, 2 );

  for (std::size_t m=0; m<2; ++m) {  // for both mesh partitions

    // Query number of nodes in mesh partition
    auto npoin = psup[m].second.size()-1;

    // Create CSR matrix based on mesh partition with psup
    tk::CSR A( /* DOF = */ 1, psup[m] );

    const auto& X = coord[m][0];
    const auto& Y = coord[m][1];
    const auto& Z = coord[m][2];

    // fill matrix for mesh partition with Laplacian
    for (std::size_t e=0; e<inpoel[m].size()/4; ++e) {
      // access node IDs
      const std::array< std::size_t, 4 >
        N{{ inpoel[m][e*4+0], inpoel[m][e*4+1],
            inpoel[m][e*4+2], inpoel[m][e*4+3] }};
      // compute element Jacobi determinant
      const std::array< tk::real, 3 >
        ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
        ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
        da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
      const auto J = tk::triple( ba, ca, da );        // J = 6V
      Assert( J > 0, "Element Jacobian non-positive" );

      // shape function derivatives, nnode*ndim [4][3]
      std::array< std::array< tk::real, 3 >, 4 > grad;
      grad[1] = tk::crossdiv( ca, da, J );
      grad[2] = tk::crossdiv( da, ba, J );
      grad[3] = tk::crossdiv( ba, ca, J );
      for (std::size_t i=0; i<3; ++i)
        grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

      for (std::size_t a=0; a<4; ++a)
        for (std::size_t k=0; k<3; ++k)
           for (std::size_t b=0; b<4; ++b)
             A(N[a],N[b]) += J/6 * grad[a][k] * grad[b][k];
    }

    // Create RHS and solution/unknown vectors for mesh partition
    std::vector< tk::real > b(npoin,1.0), x(npoin,0.0);

    // Grab a node (gid=0) as Dirichlet BC if on partition
    auto it = lid[m].find( 0 );
    if (it != end(lid[m])) { // if row with global id g exists on this partition
      auto i = it->second;
      A.dirichlet( i, gid[m], nodecommap[m] );
    }

    // Dynamically insert array element with mesh partition
    auto M = static_cast< int >( m );
    cg[M].insert( A, x, b, gid[m], lid[m], nodecommap[m] );

    // Initialize CG solve for mesh partition
    cg[M].setup(CkCallback(CkIndex_CGReceiver::initialized(nullptr), host[M]));

  }

  // Start reduction manager on CG chare array
  cg.doneInserting();
}

//! Test simple CG solve with DOF=3 of Laplacian in serial
template<> template<>
void ConjugateGradients_object::test< 3 >() {
  // This test spawns a new Charm++ chare which creates two new TUT tests. Each
  // test has similar test names but differ at the end, corresponding to a
  // stage of the completion of the multiple tests as part of the series. The
  // spawed tests create new test results, sending them back to the suite if
  // successful.
  set_test_name( "Laplacian 3DOF 1PE setup" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel {
    3, 13, 8, 14,
    12, 3, 13, 8,
    8, 3, 14, 11,
    12, 3, 8, 11,
    1, 2, 3, 13,
    6, 13, 7, 8,
    5, 9, 14, 11,
    5, 1, 3, 14,
    10, 4, 12, 11,
    2, 6, 12, 13,
    8, 7, 9, 14,
    13, 1, 7, 14,
    5, 3, 4, 11,
    6, 10, 12, 8,
    3, 2, 4, 12,
    10, 8, 9, 11,
    3, 1, 13, 14,
    13, 7, 8, 14,
    6, 12, 13, 8,
    9, 8, 14, 11,
    3, 5, 14, 11,
    4, 3, 12, 11,
    3, 2, 12, 13,
    10, 12, 8, 11 };

  // Mesh node coordinates for simple tet mesh above
  std::array< std::vector< tk::real >, 3 > coord {{
    {{ -0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0 }},
    {{ 0.5, 0.5, 0, -0.5, -0.5, 0.5, 0.5, 0, -0.5, -0.5, -0.5, 0, 0.5, 0 }},
    {{ -0.5, 0.5, 0, 0.5, -0.5, 0.5, -0.5, 0, -0.5, 0.5, 0, 0.5, 0, -0.5 }} }};

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );
  // Query number of nodes in mesh
  auto npoin = psup.second.size()-1;

  // Create CSR matrix based on mesh with psup
  tk::CSR A( /* DOF = */ 3, psup );

  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];

  // fill matrix with Laplacian
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );

    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    for (std::size_t a=0; a<4; ++a)
      for (std::size_t k=0; k<3; ++k)
         for (std::size_t b=0; b<4; ++b)
           for (std::size_t i=0; i<3; ++i)
             A(N[a],N[b],i) += J/6 * grad[a][k] * grad[b][k];
  }

  // Create RHS and solution/unknown vectors
  std::vector< tk::real > b(npoin*3,1.0), x(npoin*3,0.0);

  // Grab a node (id=0 for all components) as Dirichlet BC
  for (std::size_t i=0; i<3; ++i) A.dirichlet( 0, {}, {}, i );

  // Create CG solver (chare array with a single element)
  tk::CProxy_ConjugateGradients cg =
    tk::CProxy_ConjugateGradients::ckNew( A, x, b, {}, {}, {}, 1 );

  // Create receiver chare whose callbacks are called when a CG task is done
  CProxy_CGReceiver host =
    CProxy_CGReceiver::ckNew( "Laplacian 3DOF 1PE",
                              /* maxit = */ 1000,
                              /* tol = * */ 1.0e-15,
                              /* normb_ex = */ 6.4807406984078604,
                              /* normres_ex = */ 5.7532468516593534e-14,
                              cg, 1 );

  // Initialize CG solve
  cg[0].setup( CkCallback(CkIndex_CGReceiver::initialized(nullptr), host[0]) );
}

//! Test simple CG solve with DOF=3 of Laplacian on 2 PEs
template<> template<>
void ConjugateGradients_object::test< 4 >() {
  // This test spawns a new Charm++ chare array of two elements which each
  // creates two new TUT tests, thus resulting four newly spawned TUT tests.
  // Each test has similar test names but differ, corresponding to a stage of
  // the completion of the multiple tests as part of the series and
  // corresponding to which chare array element they complete. The spawed tests
  // create new test results, sending them back to the suite if successful.
  set_test_name( "Laplacian 3DOF 2PE setup" );

  // Initialize psup for PE0 and PE1
  std::vector<
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > > psup {
    { { 0,1,2,4,8,9,0,2,3,7,8,0,1,3,4,5,6,7,8,9,1,2,4,6,7,0,2,3,6,9,2,6,7,8,9,2,
        3,4,5,7,9,1,2,3,5,6,8,0,1,2,5,7,9,0,2,4,5,6,8 },
    { 0,5,10,19,24,29,34,40,46,52,58 } },
    { { 0,5,11,12,4,10,11,8,9,10,7,9,12,1,5,6,8,10,11,0,4,6,7,11,12,4,5,7,8,9,10,
      11,12,3,5,6,8,9,12,2,4,6,7,9,10,2,3,6,7,8,10,12,1,2,4,6,8,9,11,0,1,4,5,6,
      10,12,0,3,5,6,7,9,11 },
    { 0,3,6,9,12,18,24,32,38,44,51,58,65,72 } } };

  // Initialize node coordinates for PE0 and PE1
  std::vector< std::array< std::vector< tk::real >, 3 > > coord {{
    {{ {{ -0.5,-0.5,-0.5,-0.5,-0.5,0.5,0,0,0,0 }},
       {{ 0.5,0.5,0,-0.5,-0.5,0,-0.5,0,0.5,0 }},
       {{ -0.5,0.5,0,0.5,-0.5,0,0,0.5,0,-0.5 }} }},
    {{ {{ -0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0 }},
       {{ 0.5,0.5,-0.5,-0.5,0.5,0.5,0,-0.5,-0.5,-0.5,0,0.5,0 }},
       {{ -0.5,0.5,0.5,-0.5,0.5,-0.5,0,-0.5,0.5,0,0.5,0,-0.5 }} }} }};

  // Initialize mesh connectivity for PE0 and PE1
  std::vector< std::vector< std::size_t > > inpoel {
    { 2,8,5,9,7,2,8,5,5,2,9,6,7,2,5,6,0,1,2,8,4,0,2,9,4,2,3,6,2,1,3,7,
      2,0,8,9,2,4,9,6,3,2,7,6,2,1,7,8 },
    { 4,8,10,6,8,6,7,9,11,5,6,12,4,10,11,6,7,6,12,9,8,10,6,9,4,11,5,6,
      3,7,12,9,8,2,10,9,1,4,10,11,6,5,7,12,11,0,5,12 } };

  // Initialize global mesh node IDs for PE0 and PE1
  std::vector< std::vector< std::size_t > > gid {
    { 0,1,2,3,4,7,10,11,12,13 },
    { 0,1,3,4,5,6,7,8,9,10,11,12,13 } };

  // Initialize local->global mesh node ID map for PE0 and PE1
  std::vector< std::unordered_map< std::size_t, std::size_t > > lid{
    { {0,0},{11,7},{1,1},{12,8},{2,2},{13,9},{3,3},{4,4},{7,5},{10,6} },
    { {0,0},{11,10},{5,4},{1,1},{3,2},{4,3},{6,5},{7,6},{8,7},{9,8},
      {10,9},{12,11},{13,12} } };

  // Initialize node communication map for PE0 and PE1
  std::vector< tk::NodeCommMap > nodecommap;
  nodecommap.push_back({});
  nodecommap.back()[1].insert( { 7,4,10,12,1,13,3,0,11 } );
  nodecommap.push_back({});
  nodecommap.back()[0].insert( { 7,4,3,10,1,12,11,0,13 } );

  // Create CG solver (empty chare array, will use dynamic insertion)
  tk::CProxy_ConjugateGradients cg = tk::CProxy_ConjugateGradients::ckNew();

  // Create receiver chare array (2 elements) whose callbacks are called when a
  // CG task is done on a PE
  CProxy_CGReceiver host =
    CProxy_CGReceiver::ckNew( "Laplacian 3DOF 2PE",
                              /* maxit = */ 1000,
                              /* tol = * */ 1.0e-15,
                              /* normb_ex = */ 6.4807406984078604,
                              /* normres_ex = */ 5.7532468516593534e-14,
                              cg, 2 );

  for (std::size_t m=0; m<2; ++m) {  // for both mesh partitions

    // Query number of nodes in mesh partition
    auto npoin = psup[m].second.size()-1;

    // Create CSR matrix based on mesh partition with psup
    tk::CSR A( /* DOF = */ 3, psup[m] );

    const auto& X = coord[m][0];
    const auto& Y = coord[m][1];
    const auto& Z = coord[m][2];

    // fill matrix for mesh partition with Laplacian
    for (std::size_t e=0; e<inpoel[m].size()/4; ++e) {
      // access node IDs
      const std::array< std::size_t, 4 >
        N{{ inpoel[m][e*4+0], inpoel[m][e*4+1],
            inpoel[m][e*4+2], inpoel[m][e*4+3] }};
      // compute element Jacobi determinant
      const std::array< tk::real, 3 >
        ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
        ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
        da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
      const auto J = tk::triple( ba, ca, da );        // J = 6V
      Assert( J > 0, "Element Jacobian non-positive" );

      // shape function derivatives, nnode*ndim [4][3]
      std::array< std::array< tk::real, 3 >, 4 > grad;
      grad[1] = tk::crossdiv( ca, da, J );
      grad[2] = tk::crossdiv( da, ba, J );
      grad[3] = tk::crossdiv( ba, ca, J );
      for (std::size_t i=0; i<3; ++i)
        grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

      for (std::size_t a=0; a<4; ++a)
        for (std::size_t k=0; k<3; ++k)
           for (std::size_t b=0; b<4; ++b)
             for (std::size_t i=0; i<3; ++i)
               A(N[a],N[b],i) += J/6 * grad[a][k] * grad[b][k];
    }

    // Create RHS and solution/unknown vectors for mesh partition
    std::vector< tk::real > b(npoin*3,1.0), x(npoin*3,0.0);

    // Grab a node (gid=0 for all components) as Dirichlet BC if on partition
    auto it = lid[m].find( 0 );
    if (it != end(lid[m])) { // if row with global id g exists on this partition
      auto i = it->second;
      for (std::size_t j=0; j<3; ++j)
        A.dirichlet( i, gid[m], nodecommap[m], j );
    }

    // Dynamically insert array element with mesh partition
    auto M = static_cast< int >( m );
    cg[M].insert( A, x, b, gid[m], lid[m], nodecommap[m] );

    // Initialize CG solve for mesh partition
    cg[M].setup(CkCallback(CkIndex_CGReceiver::initialized(nullptr), host[M]));

  }

  // Start reduction manager on CG chare array
  cg.doneInserting();
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT

#include "NoWarning/cgreceiver.def.h"
