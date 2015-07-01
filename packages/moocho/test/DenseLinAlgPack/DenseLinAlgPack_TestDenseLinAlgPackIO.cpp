// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <ostream>
#include <istream>
#include <iomanip>

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"
#include "DenseLinAlgPack_IO.hpp"
#include "DenseLinAlgPack_OutFormat.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"


// To use this function you must bind _in to the file TestDenseLinAlgPackIO.in

void DenseLinAlgPack::TestingPack::TestDenseLinAlgPackIO(std::istream& _in, std::ostream& _out) {

  using std::endl;
  using std::setw;

  try {

  using InputStreamHelperPack::eat_comment_lines;
  using LinAlgPackIO::bind;
  using LinAlgPackIO::cbind;
  using LinAlgPackIO::format;
  
  _out
    << "\n*****************************************************************"
    << "\n*** Testing input/output operators and functions for DVector, ***"
    << "\n*** DVectorSlice, DMatrix and DMatrixSlice.                   ***"
    << "\n*** Must be run with TestDenseLinAlgPackIO.in as input        ***"
    << "\n*****************************************************************\n";
  
  // Creating a formating object
  _out	<< "\n  format f(_in);\n"
      << "  f.setw(15).showpoint().setprecision(6);\n";
  format f(_in);
  f.setw(15).showpoint().setprecision(6);

  // Testing The eating of comment lines
  _out	<< "\nTest Eating comment lines.\n"
      << "  eat_comment_lines(_in,'*');\n";
  eat_comment_lines(_in,'*');

  // Test inputing and outputing Vectors and VectorSlices

  _out	<< "\nTest inputing Vectors and VectorSlices (eat_comment_lines(_in,'*') called inbetween)\n"
      << "  DVector v1, v2(4), v3, v4(4);\n"
      << "\n  _in >> v1;\n";
  DVector v1, v2(4), v3, v4(4);
  _in >> v1;							// input (1) : DVector with resizing from 0
  _out	<< "\n  _out	<< cbind(f,v1);\n";
  _out	<< cbind(f,v1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> v2;\n";
  _in >> v2;							// input (2) : DVector with resizing from sized
  _out	<< "\n  _out	<< cbind(f,v2);\n";
  _out	<< cbind(f,v2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> v3;\n";
  _in >> v3;							// input (3) : DVector with resizing from 0
  _out	<< "\n  _out	<< cbind(f,v3());\n";
  _out	<< cbind(f,v3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> v4();\n";
  { DVectorSlice t = v4(); _in >> t; } // input (4) : DVectorSlice with no resizing (must be size 4)
  _out	<< "\n  _out	<< cbind(f,v4());\n";
  _out	<< cbind(f,v4());

  eat_comment_lines(_in,'*');

  _out	<< "  v1.resize(0); v3.resize(0);\n"
      << "\n  _in >> bind(f,v1);\n";
  v1.resize(0); v3.resize(0);
  _in >> bind(f,v1);				// input (5) : DVector with bound format, resize from 0
  _out	<< "\n  _out	<< cbind(f,v1);\n";
  _out	<< cbind(f,v1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v2);\n";
  _in >> bind(f,v2);				// input (6) : DVector with bound format, resize from sized
  _out	<< "\n  _out	<< cbind(f,v2);\n";
  _out	<< cbind(f,v2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v3);\n";
  _in >> bind(f,v3);				// input (7) : DVector with bound format, resize form 0
  _out	<< "\n  _out	<< cbind(f,v3());\n";
  _out	<< cbind(f,v3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v4());\n";
  { DVectorSlice t = v4(); _in >> bind(f,t); }		// input (8) : DVectorSlice with bound format, no resizing (must be size 4)
  _out	<< "\n  _out	<< cbind(f,v4());\n";
  _out	<< cbind(f,v4());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f.ignore_dim(),v1);\n";
  _in >> bind(f.ignore_dim(),v1);	// input (9) : DVector with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,v1);\n";
  _out	<< cbind(f,v1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v2);\n";
  _in >> bind(f,v2);				// input (10) : DVector with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,v2);\n";
  _out	<< cbind(f,v2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v3());\n";
  { DVectorSlice t = v3(); _in >> bind(f,t); }    // input (11) : DVectorSlice with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,v3());\n";
  _out	<< cbind(f,v3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,v4());\n";
  { DVectorSlice t = v4(); _in >> bind(f,t); }       // input (12) : DVectorSlice with bound format,  ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,v4());\n";
  _out	<< cbind(f,v4());

  // Test variations of outputing for DVector and DVectorSlice objects

  _out	<< "\nTest variations of outputing for Vectors and VectorSlices\n";

  _out	<< "\n  _out << cbind(f.ignore_dim(),v1);\n";
  _out << cbind(f.ignore_dim(),v1);

  _out	<< "\n  _out << cbind(f.no_ignore_dim(), const_cast<const DVectorSlice&>(v1()) );\n";
  { const DVectorSlice t = v1(); _out << cbind(f.no_ignore_dim(), t ); }

  _out	<< "\n  _out << cbind(f.ignore_dim().no_insert_newlines(),v1) << cbind(f,v1) << endl;\n";
  _out << cbind(f.ignore_dim().no_insert_newlines(),v1) << cbind(f,v1) << endl;

  // Test inputing and outputing DMatrix and DMatrixSlice objects

  eat_comment_lines(_in,'*');

  _out	<< "\nTest inputing and outputing DMatrix and DMatrixSlice objects (eat_comment_lines(_in,'*') called inbetween)\n"
      << "  DMatrix m1, m2(2,2), m3, m4(2,2);\n"
      << "\n  _in >> m1;\n";
  DMatrix m1, m2(2,2), m3, m4(2,2);
  _in >> m1;							// input (13) : DMatrix with resizing from 0
  _out	<< "\n  _out	<< cbind(f.no_ignore_dim().insert_newlines(),m1);\n";
  _out	<< cbind(f.no_ignore_dim().insert_newlines(),m1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> m2;\n";
  _in >> m2;							// input (14) : DMatrix with resizing from sized
  _out	<< "\n  _out	<< cbind(f,m2);\n";
  _out	<< cbind(f,m2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> m3;\n";
  _in >> m3;							// input (15) : DMatrix with resizing from 0
  _out	<< "\n  _out	<< cbind(f,m3());\n";
  _out	<< cbind(f,m3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> m4();\n";
  { DMatrixSlice t = m4(); _in >> t; }        // input (16) : DMatrixSlice with no resizing (must be size 4)
  _out	<< "\n  _out	<< cbind(f,m4());\n";
  _out	<< cbind(f,m4());

  eat_comment_lines(_in,'*');

  _out	<< "  m1.resize(0,0); m3.resize(0,0);\n"
      << "\n  _in >> bind(f,m1);\n";
  m1.resize(0,0); m3.resize(0,0);
  _in >> bind(f,m1);				// input (17) : DMatrix with bound format, resize from 0
  _out	<< "\n  _out	<< cbind(f,m1);\n";
  _out	<< cbind(f,m1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m2);\n";
  _in >> bind(f,m2);				// input (18) : DMatrix with bound format, resize from sized
  _out	<< "\n  _out	<< cbind(f,m2);\n";
  _out	<< cbind(f,m2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m3);\n";
  _in >> bind(f,m3);				// input (19) : DMatrix with bound format, resize form 0
  _out	<< "\n  _out	<< cbind(f,m3());\n";
  _out	<< cbind(f,m3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m4());\n";
  { DMatrixSlice t = m4(); _in >> bind(f,t); }   // input (20) : DMatrixSlice with bound format, no resizing (must be size 4)
  _out	<< "\n  _out	<< cbind(f,m4());\n";
  _out	<< cbind(f,m4());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f.ignore_dim(),m1);\n";
  _in >> bind(f.ignore_dim(),m1);	// input (21) : DMatrix with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,m1);\n";
  _out	<< cbind(f,m1);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m2);\n";
  _in >> bind(f,m2);				// input (22) : DMatrix with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,m2);\n";
  _out	<< cbind(f,m2);

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m3);\n";
  _in >> bind(f,m3);				// input (23) : DMatrixSlice with bound format, ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,m3());\n";
  _out	<< cbind(f,m3());

  eat_comment_lines(_in,'*');

  _out	<< "\n  _in >> bind(f,m4());\n";
  { DMatrixSlice t = m4(); _in >> bind(f,t); } // input (24) : DMatrixSlice with bound format,  ignore vector dimension
  _out	<< "\n  _out	<< cbind(f,m4());\n";
  _out	<< cbind(f,m4());

  // Test variations of outputing for DMatrix and DMatrixSlice objects

  _out	<< "\nTest variations of outputing for DMatrix and DMatrixSlice objects\n";

  _out	<< "\n  _out << cbind(f.ignore_dim(),m1);\n";
  _out << cbind(f.ignore_dim(),m1);

  _out	<< "\n  _out << cbind(f.no_ignore_dim(), const_cast<const DMatrixSlice&>(m1()) );\n";
  { const DMatrixSlice t = m1(); _out << cbind(f.no_ignore_dim(), t ); }

  _out	<< "\n  _out << 2*m1.rows() << ' ' << m1.cols() << endl << cbind(f.ignore_dim(),m1) << cbind(f,m1);\n";
  _out << 2*m1.rows() << ' ' << m1.cols() << endl << cbind(f.ignore_dim(),m1) << cbind(f,m1);

  _out	<< "\nIf you read this then no unexpected exceptions occured.\n";

  } // end try
  catch(const std::exception& excpt) {
    _out << "\nCaught a std::exception: " << excpt.what() << endl;
  }
  catch(...) {
    _out << "\nCaught and unknown exception\n";
  }
}
