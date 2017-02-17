// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_02.cpp
\brief  Unit test for the scalar multiply operations of the ArrayTools class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <Kokkos_Core.hpp>
using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Unit Test (ArrayTools)                                |\n" \
  << "|                                                                             |\n" \
  << "|     1) Array operations: scalar multiply                                    |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = Teuchos::TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 36;
#endif
  typedef ArrayTools art; 
  typedef RealSpaceTools<double> rst; 
#ifdef HAVE_INTREPID_DEBUG
  ArrayTools atools;
#endif
  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: exceptions                                                          |\n"\
  << "===============================================================================\n";

  try{

#ifdef HAVE_INTREPID_DEBUG
    FieldContainer<double> a_2_2(2, 2);
    FieldContainer<double> a_10_2(10, 2);
    FieldContainer<double> a_10_3(10, 3);
    FieldContainer<double> a_10_2_2(10, 2, 2);
    FieldContainer<double> a_10_2_3(10, 2, 3);
    FieldContainer<double> a_10_3_2(10, 3, 2);
    FieldContainer<double> a_9_2_2(9, 2, 2);

    FieldContainer<double> a_10_2_2_2(10, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2(9, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2(10, 3, 2, 2);
    FieldContainer<double> a_10_2_3_2(10, 2, 3, 2);
    FieldContainer<double> a_10_2_2_3(10, 2, 2, 3);

    FieldContainer<double> a_10_2_2_2_2(10, 2, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2_2(9, 2, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2_2(10, 3, 2, 2, 2);
    FieldContainer<double> a_10_2_3_2_2(10, 2, 3, 2, 2);
    FieldContainer<double> a_10_2_2_3_2(10, 2, 2, 3, 2);
    FieldContainer<double> a_10_2_2_2_3(10, 2, 2, 2, 3);

    FieldContainer<double> a_9_2(9, 2);
    FieldContainer<double> a_10_1(10, 1);

    FieldContainer<double> a_10_1_2(10, 1, 2);
    FieldContainer<double> a_10_1_3(10, 1, 3);

    FieldContainer<double> a_10_1_2_2(10, 1, 2, 2);

    FieldContainer<double> a_2_3_2_2(2, 3, 2, 2);
    FieldContainer<double> a_2_2_2_2(2, 2, 2, 2);
    FieldContainer<double> a_2_10(2, 10);
    FieldContainer<double> a_2(2);

    *outStream << "-> scalarMultiplyDataField:\n";
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_2_2, a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_2_2, a_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2, a_10_3, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_9_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_3_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_3_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_3_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_3, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_2, a_10_1, a_10_2_2_2_2) );
    //
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_2_2, a_10_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2, a_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2, a_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2, a_10_2, a_2_10) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_2, a_9_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_3_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_3_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_3_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_3, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataField<double>(a_10_2_2_2_2, a_10_1, a_2_2_2_2) );


    FieldContainer<double> a_2_2_2(2, 2, 2);

    *outStream << "-> scalarMultiplyDataData:\n";
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_2_2, a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_2, a_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2, a_10_3, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_9_2_2_2, a_10_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_3_2_2, a_10_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_3_2, a_10_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_3, a_10_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_10_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_10_1, a_10_2_2_2) );
    //
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_2_2, a_10_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2_2, a_2_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2, a_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_9_2, a_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_3_2_2, a_10_2, a_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_3_2, a_10_2, a_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_3, a_10_2, a_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_10_2, a_2_2_2) );
    INTREPID_TEST_COMMAND( atools.scalarMultiplyDataData<double>(a_10_2_2_2, a_10_1, a_2_2_2) );
#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

#ifdef HAVE_INTREPID_DEBUG
  if (Teuchos::TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;
#endif


  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=false (branch 1) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      Kokkos::View<double***> in_c_f_p("in_c_f_p", c, f, p);
      Kokkos::View<double****> in_c_f_p_d("in_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> in_c_f_p_d_d("in_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double***> out_c_f_p("out_c_f_p", c, f, p);
      Kokkos::View<double***> outi_c_f_p("outi_c_f_p", c, f, p);
      Kokkos::View<double****> out_c_f_p_d("out_c_f_p_d", c, f, p, d1);
      Kokkos::View<double****> outi_c_f_p_d("outi_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> out_c_f_p_d_d("out_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double*****> outi_c_f_p_d_d("outi_c_f_p_d_d", c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
  for (unsigned int i=0; i<in_c_f_p.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p.dimension(2); k++)
      in_c_f_p(i,j,k) = Teuchos::ScalarTraits<double>::random();

  for (unsigned int i=0; i<in_c_f_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d.dimension(3); l++)
        in_c_f_p_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
      
  for (unsigned int i=0; i<in_c_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d_d.dimension(3); l++)
          for (unsigned int m=0; m<in_c_f_p_d_d.dimension(4); m++)      
          in_c_f_p_d_d(i,j,k,l,m) = Teuchos::ScalarTraits<double>::random();
      
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
		}
		
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
		}		
		


      art::scalarMultiplyDataField<double>(out_c_f_p, data_c_p, in_c_f_p);
      art::scalarMultiplyDataField<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::vectorNorm(outi_c_f_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::vectorNorm(outi_c_f_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
 for (unsigned int i=0; i<in_c_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d_d.dimension(3); l++)
          for (unsigned int m=0; m<in_c_f_p_d_d.dimension(4); m++)      
          in_c_f_p_d_d(i,j,k,l,m) = 1.0;
             
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
		}
		
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
		}	

      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_p(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_1(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=false (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      Kokkos::View<double**> in_f_p("in_f_p", f, p);
      Kokkos::View<double***> in_f_p_d("in_f_p_d", f, p, d1);
      Kokkos::View<double****> in_f_p_d_d("in_f_p_d_d", f, p, d1, d2);
      Kokkos::View<double***> in_c_f_p("in_c_f_p", c, f, p);
      Kokkos::View<double****> in_c_f_p_d("in_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> in_c_f_p_d_d("in_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> data_c_p_one("data_c_p_one", c, p);
      Kokkos::View<double**> data_c_1_one("data_c_1_one", c, 1);
      Kokkos::View<double***> out_c_f_p("out_c_f_p", c, f, p);
      Kokkos::View<double***> outi_c_f_p("outi_c_f_p", c, f, p);
      Kokkos::View<double****> out_c_f_p_d("out_c_f_p_d", c, f, p, d1);
      Kokkos::View<double****> outi_c_f_p_d("outi_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> out_c_f_p_d_d("out_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double*****> outi_c_f_p_d_d("outi_c_f_p_d_d", c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
  for (unsigned int i=0; i<in_f_p.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p.dimension(1); j++)
      in_f_p(i,j) = Teuchos::ScalarTraits<double>::random();
      
  for (unsigned int i=0; i<in_f_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d.dimension(2); k++)
      in_f_p_d(i,j,k) = Teuchos::ScalarTraits<double>::random();      
  
  for (unsigned int i=0; i<in_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_f_p_d_d.dimension(3); l++)      
          in_f_p_d_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();

  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
		data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
        data_c_p_one(i,j) = 1.0;	
	}

  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
		data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
        data_c_1_one(i,j) = 1.0;
		
	}

      art::scalarMultiplyDataField<double>(out_c_f_p, data_c_p, in_f_p);
      art::scalarMultiplyDataField<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      art::scalarMultiplyDataField<double>(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::vectorNorm(outi_c_f_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d, data_c_p, in_f_p_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      art::scalarMultiplyDataField<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::vectorNorm(outi_c_f_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants

  for (unsigned int i=0; i<in_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_f_p_d_d.dimension(3); l++)      
          in_f_p_d_d(i,j,k,l) = 1.0;      
      
      
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++)
      data_c_p(i,j) = 5.0;
      
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++)
      data_c_1(i,j) = 5.0;      

      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_p(0,0)*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_f_p_d_d(0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_1(0,0)*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_1(0,0)*in_f_p_d_d(0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=true, i.e. division (branch 1) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      Kokkos::View<double***> in_c_f_p("in_c_f_p", c, f, p);
      Kokkos::View<double****> in_c_f_p_d("in_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> in_c_f_p_d_d("in_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double***> out_c_f_p("out_c_f_p", c, f, p);
      Kokkos::View<double***> outi_c_f_p("outi_c_f_p", c, f, p);
      Kokkos::View<double****> out_c_f_p_d("out_c_f_p_d", c, f, p, d1);
      Kokkos::View<double****> outi_c_f_p_d("outi_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> out_c_f_p_d_d("out_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double*****> outi_c_f_p_d_d("outi_c_f_p_d_d", c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
   for (unsigned int i=0; i<in_c_f_p.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p.dimension(2); k++)
      in_c_f_p(i,j,k) = Teuchos::ScalarTraits<double>::random(); 
      
  for (unsigned int i=0; i<in_c_f_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d.dimension(3); l++)      
          in_c_f_p_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
      
 for (unsigned int i=0; i<in_c_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d_d.dimension(3); l++)
          for (unsigned int m=0; m<in_c_f_p_d_d.dimension(4); m++)  
           in_c_f_p_d_d(i,j,k,l,m) = Teuchos::ScalarTraits<double>::random();

  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);

}
 
   for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);

}
   

      art::scalarMultiplyDataField<double>(out_c_f_p, data_c_p, in_c_f_p, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p, datainv_c_p, out_c_f_p, true);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::vectorNorm(outi_c_f_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d, data_c_p, in_c_f_p_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d, true);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::vectorNorm(outi_c_f_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d, true);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d, true);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
 for (unsigned int i=0; i<in_c_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_f_p_d_d.dimension(3); l++)
          for (unsigned int m=0; m<in_c_f_p_d_d.dimension(4); m++)  
           in_c_f_p_d_d(i,j,k,l,m) = 1.0;      
      
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;


}  

   for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;

}

      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d, true);
       if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2)/rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2)/rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d,  NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      
      } // end scope
     { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=true, i.e. division (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      Kokkos::View<double**> in_f_p("in_f_p", f, p);
      Kokkos::View<double***> in_f_p_d("in_f_p_d", f, p, d1);
      Kokkos::View<double****> in_f_p_d_d("in_f_p_d_d", f, p, d1, d2);
      Kokkos::View<double***> in_c_f_p("in_c_f_p", c, f, p);
      Kokkos::View<double****> in_c_f_p_d("in_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> in_c_f_p_d_d("in_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> data_c_p_one("data_c_p_one", c, p);
      Kokkos::View<double**> data_c_1_one("data_c_1_one", c, 1);
      Kokkos::View<double***> out_c_f_p("out_c_f_p", c, f, p);
      Kokkos::View<double***> outi_c_f_p("outi_c_f_p", c, f, p);
      Kokkos::View<double****> out_c_f_p_d("out_c_f_p_d", c, f, p, d1);
      Kokkos::View<double****> outi_c_f_p_d("outi_c_f_p_d", c, f, p, d1);
      Kokkos::View<double*****> out_c_f_p_d_d("out_c_f_p_d_d", c, f, p, d1, d2);
      Kokkos::View<double*****> outi_c_f_p_d_d("outi_c_f_p_d_d", c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
  for (unsigned int i=0; i<in_f_p.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p.dimension(1); j++)
        in_f_p(i,j) = Teuchos::ScalarTraits<double>::random();
           

  for (unsigned int i=0; i<in_f_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d.dimension(2); k++)
       in_f_p_d(i,j,k) = Teuchos::ScalarTraits<double>::random();
                  
  for (unsigned int i=0; i<in_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_f_p_d_d.dimension(3); l++)                   
        in_f_p_d_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
     
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
        data_c_p_one(i,j) = 1.0;
      	 }    
     
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
        data_c_1_one(i,j) = 1.0;
      	 }                     

      art::scalarMultiplyDataField<double>(out_c_f_p, data_c_p, in_f_p, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p, datainv_c_p, out_c_f_p, true);
      art::scalarMultiplyDataField<double>(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::vectorNorm(outi_c_f_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d, data_c_p, in_f_p_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d, true);
      art::scalarMultiplyDataField<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::vectorNorm(outi_c_f_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d, true);
      art::scalarMultiplyDataField<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      
   for (unsigned int i=0; i<in_f_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_f_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_f_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_f_p_d_d.dimension(3); l++)                   
        in_f_p_d_d(i,j,k,l) = 1.0;
        
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
      	 }            

  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
      	 }        

      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2)/rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2)/rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_1(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope
 
      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=false (branch 1) ************\n";

      int c=5, p=9, d1=7, d2=13;

      Kokkos::View<double**> in_c_p("in_c_p", c, p);
      Kokkos::View<double***> in_c_p_d("in_c_p_d", c, p, d1);
      Kokkos::View<double****> in_c_p_d_d("in_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> out_c_p("out_c_p", c, p);
      Kokkos::View<double**> outi_c_p("outi_c_p", c, p);
      Kokkos::View<double***> out_c_p_d("out_c_p_d", c, p, d1);
      Kokkos::View<double***> outi_c_p_d("outi_c_p_d", c, p, d1);
      Kokkos::View<double****> out_c_p_d_d("out_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double****> outi_c_p_d_d("outi_c_p_d_d", c, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
  for (unsigned int i=0; i<in_c_p.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p.dimension(1); j++)
        in_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
      
  for (unsigned int i=0; i<in_c_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d.dimension(2); k++)
       in_c_p_d(i,j,k) = Teuchos::ScalarTraits<double>::random();      

  for (unsigned int i=0; i<in_c_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_p_d_d.dimension(3); l++)                   
        in_c_p_d_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
        
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
         datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
      	 }           
      	  
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
      	 }         

      art::scalarMultiplyDataData<double>(out_c_p, data_c_p, in_c_p);
      art::scalarMultiplyDataData<double>(outi_c_p, datainv_c_p, out_c_p);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::vectorNorm(outi_c_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d, data_c_p, in_c_p_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d, datainv_c_p, out_c_p_d);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::vectorNorm(outi_c_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_c_p_d_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_c_p_d_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
  for (unsigned int i=0; i<in_c_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_p_d_d.dimension(3); l++)                   
        in_c_p_d_d(i,j,k,l) = 1.0;
 
      
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
      	 }           

  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
      	 }     
      	
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_c_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_c_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_1(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

     { // start scope

 //     *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=false (branch 2) ************\n";
 

      int c=5, p=9, d1=7, d2=13;

      Kokkos::View<double*> in_p("in_p", p);
      Kokkos::View<double**> in_p_d("in_p_d", p, d1);
      Kokkos::View<double***> in_p_d_d("in_p_d_d", p, d1, d2);
      Kokkos::View<double**> in_c_p("in_c_p", c, p);
      Kokkos::View<double***> in_c_p_d("in_c_p_d", c, p, d1);
      Kokkos::View<double****> in_c_p_d_d("in_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> data_c_p_one("data_c_p_one", c, p);
      Kokkos::View<double**> data_c_1_one("data_c_1_one", c, 1);
      Kokkos::View<double**> out_c_p("out_c_p", c, p);
      Kokkos::View<double**> outi_c_p("outi_c_p", c, p);
      Kokkos::View<double***> out_c_p_d("out_c_p_d", c, p, d1);
      Kokkos::View<double***> outi_c_p_d("outi_c_p_d", c, p, d1);
      Kokkos::View<double****> out_c_p_d_d("out_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double****> outi_c_p_d_d("outi_c_p_d_d", c, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
   for (unsigned int i=0; i<in_p.dimension(0); i++) {
        in_p(i) = Teuchos::ScalarTraits<double>::random();
      }    

  for (unsigned int i=0; i<in_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d.dimension(1); j++){
        in_p_d(i,j) = Teuchos::ScalarTraits<double>::random();
      	 }   
      	      
   for (unsigned int i=0; i<in_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_p_d_d.dimension(2); k++)
         in_p_d_d(i,j,k) = Teuchos::ScalarTraits<double>::random();

  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
        data_c_p_one(i,j) = 1.0;
      	 }   
 
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
        data_c_1_one(i,j) = 1.0;
      	 }    
 


      art::scalarMultiplyDataData<double>(out_c_p, data_c_p, in_p);
      art::scalarMultiplyDataData<double>(outi_c_p, datainv_c_p, out_c_p);
      art::scalarMultiplyDataData<double>(in_c_p, data_c_p_one, in_p);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::vectorNorm(outi_c_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d, data_c_p, in_p_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d, datainv_c_p, out_c_p_d);
      art::scalarMultiplyDataData<double>(in_c_p_d, data_c_p_one, in_p_d);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::vectorNorm(outi_c_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_p_d_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
      art::scalarMultiplyDataData<double>(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_p_d_d);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
      art::scalarMultiplyDataData<double>(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
   for (unsigned int i=0; i<in_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_p_d_d.dimension(2); k++)
         in_p_d_d(i,j,k) = 1.0;

  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
      	 }   
 
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
      	 }    
      	      
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_p(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_p_d_d);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_1(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_1(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope
      
           { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=true, i.e. division (branch 1) ************\n";

      int c=5, p=9, d1=7, d2=13;

      Kokkos::View<double**> in_c_p("in_c_p", c, p);
      Kokkos::View<double***> in_c_p_d("in_c_p_d", c, p, d1);
      Kokkos::View<double****> in_c_p_d_d("in_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> out_c_p("out_c_p", c, p);
      Kokkos::View<double**> outi_c_p("outi_c_p", c, p);
      Kokkos::View<double***> out_c_p_d("out_c_p_d", c, p, d1);
      Kokkos::View<double***> outi_c_p_d("outi_c_p_d", c, p, d1);
      Kokkos::View<double****> out_c_p_d_d("out_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double****> outi_c_p_d_d("outi_c_p_d_d", c, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      
  for (unsigned int i=0; i<in_c_p.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p.dimension(1); j++){
        in_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
      	 }   

   for (unsigned int i=0; i<in_c_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d.dimension(2); k++)
         in_c_p_d(i,j,k) = Teuchos::ScalarTraits<double>::random();

   for (unsigned int i=0; i<in_c_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_p_d_d.dimension(3); l++)    
         in_c_p_d_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
         
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
      	 }   
 
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
      	 }    
      	               

      art::scalarMultiplyDataData<double>(out_c_p, data_c_p, in_c_p, true);
      art::scalarMultiplyDataData<double>(outi_c_p, datainv_c_p, out_c_p, true);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::vectorNorm(outi_c_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d, data_c_p, in_c_p_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d, datainv_c_p, out_c_p_d, true);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::vectorNorm(outi_c_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_c_p_d_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_p, out_c_p_d_d, true);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_c_p_d_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_1, out_c_p_d_d, true);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
  for (unsigned int i=0; i<in_c_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_c_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_c_p_d_d.dimension(2); k++)
        for (unsigned int l=0; l<in_c_p_d_d.dimension(3); l++)    
         in_c_p_d_d(i,j,k,l) = 1.0;
         
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
      	 }   
 
  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
      	 }    
               
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_c_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2)/rst::vectorNorm(out_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_c_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2)/rst::vectorNorm(out_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope
      
           { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=true, i.e. division (branch 2) ************\n";

      int c=5, p=9, d1=7, d2=13;

      Kokkos::View<double*> in_p("in_p", p);
      Kokkos::View<double**> in_p_d("in_p_d", p, d1);
      Kokkos::View<double***> in_p_d_d("in_p_d_d", p, d1, d2);
      Kokkos::View<double**> in_c_p("in_c_p", c, p);
      Kokkos::View<double***> in_c_p_d("in_c_p_d", c, p, d1);
      Kokkos::View<double****> in_c_p_d_d("in_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double**> data_c_p("data_c_p", c, p);
      Kokkos::View<double**> datainv_c_p("datainv_c_p", c, p);
      Kokkos::View<double**> data_c_1("data_c_1", c, 1);
      Kokkos::View<double**> datainv_c_1("datainv_c_1", c, 1);
      Kokkos::View<double**> data_c_p_one("data_c_p_one", c, p);
      Kokkos::View<double**> data_c_1_one("data_c_1_one", c, 1);
      Kokkos::View<double**> out_c_p("out_c_p", c, p);
      Kokkos::View<double**> outi_c_p("outi_c_p", c, p);
      Kokkos::View<double***> out_c_p_d("out_c_p_d", c, p, d1);
      Kokkos::View<double***> outi_c_p_d("outi_c_p_d", c, p, d1);
      Kokkos::View<double****> out_c_p_d_d("out_c_p_d_d", c, p, d1, d2);
      Kokkos::View<double****> outi_c_p_d_d("outi_c_p_d_d", c, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
  for (unsigned int i=0; i<in_p.dimension(0); i++)
       in_p(i) = Teuchos::ScalarTraits<double>::random();
      	   
  for (unsigned int i=0; i<in_p_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d.dimension(1); j++){
        in_p_d(i,j) = Teuchos::ScalarTraits<double>::random();
      	 }         

  for (unsigned int i=0; i<in_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_p_d_d.dimension(2); k++)
         in_p_d_d(i,j,k) = Teuchos::ScalarTraits<double>::random();

  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_p(i,j) = 1.0 / data_c_p(i,j);
        data_c_p_one(i,j) = 1.0;
      }

  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = Teuchos::ScalarTraits<double>::random();
        datainv_c_1(i,j) = 1.0 / data_c_1(i,j);
        data_c_1_one(i,j) = 1.0;
      }

      art::scalarMultiplyDataData<double>(out_c_p, data_c_p, in_p, true);
      art::scalarMultiplyDataData<double>(outi_c_p, datainv_c_p, out_c_p, true);
      art::scalarMultiplyDataData<double>(in_c_p, data_c_p_one, in_p);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::vectorNorm(outi_c_p, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d, data_c_p, in_p_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d, datainv_c_p, out_c_p_d, true);
      art::scalarMultiplyDataData<double>(in_c_p_d, data_c_p_one, in_p_d);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::vectorNorm(outi_c_p_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_p_d_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_p, out_c_p_d_d, true);
      art::scalarMultiplyDataData<double>(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_p_d_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d_d, datainv_c_1, out_c_p_d_d, true);
      art::scalarMultiplyDataData<double>(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::vectorNorm(outi_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
  for (unsigned int i=0; i<in_p_d_d.dimension(0); i++)
    for (unsigned int j=0; j<in_p_d_d.dimension(1); j++)
      for (unsigned int k=0; k<in_p_d_d.dimension(2); k++){
        in_p_d_d(i,j,k) = 1.0;
      }
      
  for (unsigned int i=0; i<data_c_p.dimension(0); i++)
    for (unsigned int j=0; j<data_c_p.dimension(1); j++){
        data_c_p(i,j) = 5.0;
      }      

  for (unsigned int i=0; i<data_c_1.dimension(0); i++)
    for (unsigned int j=0; j<data_c_1.dimension(1); j++){
        data_c_1(i,j) = 5.0;
      }
      
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_p, in_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2)/rst::vectorNorm(out_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData<double>(out_c_p_d_d, data_c_1, in_p_d_d, true);
      if (std::abs(rst::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2)/rst::vectorNorm(out_c_p_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_1(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope
      /******************************************/
      *outStream << "\n";
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  Kokkos::finalize();
  return errorFlag;
}
