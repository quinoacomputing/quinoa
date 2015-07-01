/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"
#include "CEpetra_CrsGraph.h"
#include "CEpetra_Export.h"
#include "CEpetra_Import.h"
#include "Epetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID );
 **********************************************************************/

/**********************************************************************
void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );
 **********************************************************************/

/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

