/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include "Epetra_Import.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"

#include <algorithm>
#include <vector>

//==============================================================================
// Epetra_Import constructor for a Epetra_BlockMap object
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(targetMap),
    SourceMap_(sourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{

  int i;
  
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = sourceMap.NumMyElements();
  int NumTargetIDs = targetMap.NumMyElements();
  
  int *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int[NumTargetIDs];
    targetMap.MyGlobalElements(TargetGIDs);
  }
  
  int * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int[NumSourceIDs];
    sourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);
  
  
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  
  // Find count of Target IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (sourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
  
  
  
  // Define remote and permutation lists
  
  int * RemoteGIDs=0;
  RemoteLIDs_ = 0;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (sourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = sourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=targetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }

  if( NumRemoteIDs_>0 && !sourceMap.DistributedGlobal() )
    ReportError("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)", 1);
  
  // Test for distributed cases
  
  int * RemotePIDs = 0;

  if (sourceMap.DistributedGlobal()) {
    
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];
    int ierr = sourceMap.RemoteIDList(NumRemoteIDs_, RemoteGIDs, RemotePIDs, 0); // Get remote PIDs
    if (ierr) throw ReportError("Error in sourceMap.RemoteIDList call", ierr);

    //Get rid of IDs that don't exist in SourceMap
    if(NumRemoteIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumRemoteIDs_; ++i )
        if( RemotePIDs[i] == -1 ) ++cnt;
      if( cnt ) {
        if( NumRemoteIDs_-cnt ) {
          int * NewRemoteGIDs = new int[NumRemoteIDs_-cnt];
          int * NewRemotePIDs = new int[NumRemoteIDs_-cnt];
          int * NewRemoteLIDs = new int[NumRemoteIDs_-cnt];
          cnt = 0;
          for( i = 0; i < NumRemoteIDs_; ++i )
            if( RemotePIDs[i] != -1 ) {
              NewRemoteGIDs[cnt] = RemoteGIDs[i];
              NewRemotePIDs[cnt] = RemotePIDs[i];
              NewRemoteLIDs[cnt] = targetMap.LID(RemoteGIDs[i]);
              ++cnt;
            }
          NumRemoteIDs_ = cnt;
          delete [] RemoteGIDs;
          delete [] RemotePIDs;
          delete [] RemoteLIDs_;
          RemoteGIDs = NewRemoteGIDs;
          RemotePIDs = NewRemotePIDs;
          RemoteLIDs_ = NewRemoteLIDs;
          ReportError("Warning in Epetra_Import: Target IDs not found in Source Map (Do you want to import to subset of Target Map?)", 1);
        }
        else { //valid RemoteIDs empty
          NumRemoteIDs_ = 0;
          delete [] RemoteGIDs;
          RemoteGIDs = 0;
          delete [] RemotePIDs;
          RemotePIDs = 0;
        }
      }
    }

    //Sort Remote IDs by processor so DoReverses will work
    Epetra_Util util;
    int * tmpPtr[2];
    tmpPtr[0] = RemoteLIDs_, tmpPtr[1] = RemoteGIDs;
    util.Sort(true,NumRemoteIDs_,RemotePIDs,0,0,2,tmpPtr);

    Distor_ = sourceMap.Comm().CreateDistributor();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    bool Deterministic = true;
    ierr = Distor_->CreateFromRecvs( NumRemoteIDs_, RemoteGIDs, RemotePIDs,
                       Deterministic, NumExportIDs_, ExportLIDs_, ExportPIDs_ );
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);

    // Export IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumExportIDs_; i++) {
      if (ExportPIDs_[i] < 0) throw ReportError("targetMap requested a GID that is not in the sourceMap.", -1);
      ExportLIDs_[i] = sourceMap.LID(ExportLIDs_[i]);
      NumSend_ += sourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
    }
  }

  if( NumRemoteIDs_>0 ) delete [] RemoteGIDs;
  if( NumRemoteIDs_>0 ) delete [] RemotePIDs;

  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Epetra_Import copy constructor 
Epetra_Import::Epetra_Import(const Epetra_Import & Importer)
  : Epetra_Object(Importer),
    TargetMap_(Importer.TargetMap_),
    SourceMap_(Importer.SourceMap_),
    NumSameIDs_(Importer.NumSameIDs_),
    NumPermuteIDs_(Importer.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Importer.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Importer.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Importer.NumSend_),
    NumRecv_(Importer.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Importer.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Importer.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Importer.RemoteLIDs_[i];
  }

  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Importer.ExportLIDs_[i];
      ExportPIDs_[i] = Importer.ExportPIDs_[i];
    }
  }

  if (Importer.Distor_!=0) Distor_ = Importer.Distor_->Clone();

}

//==============================================================================
// Epetra_Import destructor 
Epetra_Import::~Epetra_Import()
{
  if( Distor_ != 0 ) delete Distor_;
  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by Distor_
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}
//=============================================================================
void Epetra_Import::Print(ostream & os) const
{
  // mfh 14 Dec 2011: The implementation of Print() I found here
  // previously didn't print much at all, and it included a message
  // saying that it wasn't finished ("Epetra_Import::Print needs
  // attention!!!").  What you see below is a port of
  // Tpetra::Import::print, which does have a full implementation.
  // This should allow a side-by-side comparison of Epetra_Import with
  // Tpetra::Import.

  // If true, then copy the array data and sort it before printing.
  // Otherwise, leave the data in its original order.  
  //
  // NOTE: Do NOT sort the arrays in place!  Only sort in the copy.
  // Epetra depends on the order being preserved, and some arrays'
  // orders are coupled.
  const bool sortIDs = false;

  const Epetra_Comm& comm = SourceMap_.Comm();
  const int myRank = comm.MyPID();
  const int numProcs = comm.NumProc();
  
  if (myRank == 0) {
    os << "Import Data Members:" << endl;
  }
  // We don't need a barrier before this for loop, because Proc 0 is
  // the first one to do anything in the for loop anyway.
  for (int p = 0; p < numProcs; ++p) {
    if (myRank == p) {
      os << "Image ID       : " << myRank << endl;

      os << "permuteFromLIDs:";
      if (PermuteFromLIDs_ == NULL) {
	os << " NULL";
      } else {
	std::vector<int> permuteFromLIDs (NumPermuteIDs_);
	std::copy (PermuteFromLIDs_, PermuteFromLIDs_ + NumPermuteIDs_, 
		   permuteFromLIDs.begin());
	if (sortIDs) {
	  std::sort (permuteFromLIDs.begin(), permuteFromLIDs.end());
	}
	os << " {";
	for (int i = 0; i < NumPermuteIDs_; ++i) {
	  os << permuteFromLIDs[i];
	  if (i < NumPermuteIDs_ - 1) {
	    os << ", ";
	  }
	}
	os << "}";
      }
      os << endl;

      os << "permuteToLIDs  :";
      if (PermuteToLIDs_ == NULL) {
	os << " NULL";
      } else {
	std::vector<int> permuteToLIDs (NumPermuteIDs_);
	std::copy (PermuteToLIDs_, PermuteToLIDs_ + NumPermuteIDs_, 
		   permuteToLIDs.begin());
	if (sortIDs) {
	  std::sort (permuteToLIDs.begin(), permuteToLIDs.end());
	}
	os << " {";
	for (int i = 0; i < NumPermuteIDs_; ++i) {
	  os << permuteToLIDs[i];
	  if (i < NumPermuteIDs_ - 1) {
	    os << ", ";
	  }
	}
	os << "}";
      }
      os << endl;

      os << "remoteLIDs     :";
      if (RemoteLIDs_ == NULL) {
	os << " NULL";
      } else {
	std::vector<int> remoteLIDs (NumRemoteIDs_);
	std::copy (RemoteLIDs_, RemoteLIDs_ + NumRemoteIDs_, 
		   remoteLIDs.begin());
	if (sortIDs) {
	  std::sort (remoteLIDs.begin(), remoteLIDs.end());
	}
	os << " {";
	for (int i = 0; i < NumRemoteIDs_; ++i) {
	  os << remoteLIDs[i];
	  if (i < NumRemoteIDs_ - 1) {
	    os << ", ";
	  }
	}
	os << "}";
      }
      os << endl;

      // If sorting for output, the export LIDs and export PIDs have
      // to be sorted together.  We can use Epetra_Util::Sort, using
      // the PIDs as the keys to match Tpetra::Import.
      std::vector<int> exportLIDs (NumExportIDs_);
      std::vector<int> exportPIDs (NumExportIDs_);
      if (ExportLIDs_ != NULL) {
	std::copy (ExportLIDs_, ExportLIDs_ + NumExportIDs_, exportLIDs.begin());
	std::copy (ExportPIDs_, ExportPIDs_ + NumExportIDs_, exportPIDs.begin());

	if (sortIDs && NumExportIDs_ > 0) {
	  int* intCompanions[1]; // Input for Epetra_Util::Sort().
	  intCompanions[0] = &exportLIDs[0];
	  Epetra_Util::Sort (true, NumExportIDs_, &exportPIDs[0], 
			     0, (double**) NULL, 1, intCompanions);
	}
      }

      os << "exportLIDs     :";
      if (ExportLIDs_ == NULL) {
	os << " NULL";
      } else {
	os << " {";
	for (int i = 0; i < NumExportIDs_; ++i) {
	  os << exportLIDs[i];
	  if (i < NumExportIDs_ - 1) {
	    os << ", ";
	  }
	}
	os << "}";
      }
      os << endl;

      os << "exportImageIDs :";
      if (ExportPIDs_ == NULL) {
	os << " NULL";
      } else {
	os << " {";
	for (int i = 0; i < NumExportIDs_; ++i) {
	  os << exportPIDs[i];
	  if (i < NumExportIDs_ - 1) {
	    os << ", ";
	  }
	}
	os << "}";
      }
      os << endl;

      os << "numSameIDs     : " << NumSameIDs_ << endl;
      os << "numPermuteIDs  : " << NumPermuteIDs_ << endl;
      os << "numRemoteIDs   : " << NumRemoteIDs_ << endl;
      os << "numExportIDs   : " << NumExportIDs_ << endl;

      // Epetra keeps NumSend_ and NumRecv_, whereas in Tpetra, these
      // are stored in the Distributor object.  This is why we print
      // them here.
      os << "Number of sends: " << NumSend_ << endl;
      os << "Number of recvs: " << NumRecv_ << endl;
    } // if my rank is p

    // A few global barriers give I/O a chance to complete.
    comm.Barrier();
    comm.Barrier();
    comm.Barrier();
  } // for each rank p

  const bool printMaps = false;
  if (printMaps) {
    // The original implementation printed the Maps first.  We moved
    // printing the Maps to the end, for easy comparison with the
    // output of Tpetra::Import::print().
    if (myRank == 0) {
      os << endl << endl << "Source Map:" << endl << std::flush;
    }
    comm.Barrier();
    SourceMap_.Print(os);
    comm.Barrier();
  
    if (myRank == 0) {
      os << endl << endl << "Target Map:" << endl << std::flush;
    }
    comm.Barrier();
    TargetMap_.Print(os);
    comm.Barrier();
  }

  if (myRank == 0) {
    os << endl << endl << "Distributor:" << endl << std::flush;
  }
  comm.Barrier();
  if (Distor_ == NULL) {
    if (myRank == 0) {
      os << " is NULL." << endl;
    }
  } else {
    Distor_->Print(os); // Printing the Distributor is itself distributed.
  }
  comm.Barrier();
}

