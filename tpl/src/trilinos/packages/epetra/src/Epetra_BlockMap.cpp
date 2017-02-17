
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

#include "Epetra_BlockMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_HashTable.h"

// Use the new LID hash table approach by default
#define EPETRA_BLOCKMAP_NEW_LID

//==============================================================================
// Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobal_Elements, int Element_Size, int Index_Base, const Epetra_Comm& comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(0)
{
  
  // Each processor gets roughly numGlobalPoints/p points
  // This routine automatically defines a linear partitioning of a
  // map with numGlobalPoints across the processors
  // specified in the given Epetra_Comm
  
  if (NumGlobal_Elements < 0) 
    throw ReportError("NumGlobal_Elements = " + toString(NumGlobal_Elements) + ".  Should be >= 0.", -1);
  if (Element_Size <= 0) 
    throw ReportError("ElementSize = " + toString(Element_Size) + ".  Should be > 0.", -2);
  
  BlockMapData_ = new Epetra_BlockMapData(NumGlobal_Elements, Element_Size, Index_Base, comm);
  int NumProc = comm.NumProc();
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = true;

  int MyPID = comm.MyPID();
  BlockMapData_->NumMyElements_ = BlockMapData_->NumGlobalElements_ / NumProc;
  int remainder = BlockMapData_->NumGlobalElements_ % NumProc;
  int start_index = MyPID * (BlockMapData_->NumMyElements_ + 1);

  if (MyPID < remainder) 
    BlockMapData_->NumMyElements_++;
  else 
    start_index -= (MyPID - remainder);

  BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
  BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;

  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;

  BlockMapData_->MinAllGID_ = BlockMapData_->IndexBase_;
  BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
  BlockMapData_->MinMyGID_ = start_index + BlockMapData_->IndexBase_;
  BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(BlockMapData_->NumGlobalElements_, BlockMapData_->NumMyElements_);

  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined linear distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobal_Elements, int NumMy_Elements, 
				 int Element_Size, int Index_Base, const Epetra_Comm& comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(0)
{
  if (NumGlobal_Elements < -1) 
    throw ReportError("NumGlobal_Elements = " + toString(NumGlobal_Elements) + ".  Should be >= -1.", -1);
  if (NumMy_Elements < 0) 
    throw ReportError("NumMy_Elements = " + toString(NumMy_Elements) + ".  Should be >= 0.", -2);
  if (Element_Size <= 0) 
    throw ReportError("ElementSize = " + toString(Element_Size) + ". Should be > 0.", -3);

  BlockMapData_ = new Epetra_BlockMapData(NumGlobal_Elements, Element_Size, Index_Base, comm);
  BlockMapData_->NumMyElements_ = NumMy_Elements;
  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = true;

  // Each processor gets NumMyElements points
  

  // Get processor information

  int NumProc = comm.NumProc();

  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobal_Elements, NumMy_Elements);

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc==1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobal_Elements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->	MinAllGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  }
  else if (NumProc > 1) {
    // Sum up all local element counts to get global count
    BlockMapData_->Comm_->SumAll(&BlockMapData_->NumMyElements_, &BlockMapData_->NumGlobalElements_, 1);
    
    CheckValidNGE(NumGlobal_Elements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
    
    // Use the ScanSum function to compute a prefix sum of the number of points
    BlockMapData_->Comm_->ScanSum(&BlockMapData_->NumMyElements_, &BlockMapData_->MaxMyGID_, 1);
    
    int start_index = BlockMapData_->MaxMyGID_ - BlockMapData_->NumMyElements_;
    BlockMapData_->MinMyGID_ = start_index + BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  

  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobal_Elements, int NumMy_Elements,
                                 const int * myGlobalElements, 
				 int Element_Size, int indexBase,
                                 const Epetra_Comm& comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(0)
{
  int i;
  // Each processor gets NumMyElements points

  if (NumGlobal_Elements < -1) 
    throw ReportError("NumGlobal_Elements = " + toString(NumGlobal_Elements) + ".  Should be >= -1.", -1);
  if (NumMy_Elements < 0) 
    throw ReportError("NumMy_Elements = " + toString(NumMy_Elements) + ".  Should be >= 0.", -2);
  if (Element_Size <= 0) 
    throw ReportError("ElementSize = " + toString(Element_Size) + ". Should be > 0.", -3);

  // Allocate storage for global index list information

  BlockMapData_ = new Epetra_BlockMapData(NumGlobal_Elements, Element_Size, indexBase, comm);
  if (NumMy_Elements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(NumMy_Elements);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
  }

  BlockMapData_->NumMyElements_ = NumMy_Elements;
  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = false;
  // Get processor information

  int NumProc = comm.NumProc();
  if (NumMy_Elements > 0) {
    // Compute min/max GID on this processor
    BlockMapData_->MinMyGID_ = myGlobalElements[0];
    BlockMapData_->MaxMyGID_ = myGlobalElements[0];
    for (i = 0; i < NumMy_Elements; i++) {
      BlockMapData_->MyGlobalElements_[i] = myGlobalElements[i];
      BlockMapData_->MinMyGID_ = EPETRA_MIN(BlockMapData_->MinMyGID_,myGlobalElements[i]);
      BlockMapData_->MaxMyGID_ = EPETRA_MAX(BlockMapData_->MaxMyGID_,myGlobalElements[i]);
    }
  }
  else {
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->IndexBase_ - 1;
  }
	
  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobal_Elements, NumMy_Elements);

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc==1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobal_Elements);
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->MinMyGID_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MaxMyGID_;
  }
  else if (NumProc > 1) {
    // Sum up all local element counts to get global count
    BlockMapData_->Comm_->SumAll(&BlockMapData_->NumMyElements_, &BlockMapData_->NumGlobalElements_, 1);
    CheckValidNGE(NumGlobal_Elements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    // Use the Allreduce function to find min/max GID 
    int *tmp_send = new int[2];
    int *tmp_recv = new int[2];
    tmp_send[0] = - BlockMapData_->MinMyGID_; // Negative sign lets us do one reduction
    tmp_send[1] =   BlockMapData_->MaxMyGID_;
    BlockMapData_->Comm_->MaxAll(tmp_send, tmp_recv, 2);
    BlockMapData_->MinAllGID_ = - tmp_recv[0];
    BlockMapData_->MaxAllGID_ =   tmp_recv[1];
    delete [] tmp_send;
    delete [] tmp_recv;
    if (BlockMapData_->MinAllGID_ < BlockMapData_->IndexBase_)
      throw ReportError("Minimum global element index = " + toString(BlockMapData_->MinAllGID_) + 
			" is less than index base = " + toString(BlockMapData_->IndexBase_) +".", -5);
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  

  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of variable size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobal_Elements, int NumMy_Elements,
                                 const int * myGlobalElements, 
				 const int *elementSizeList, int indexBase,
                                 const Epetra_Comm& comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(0)
{

  int i;
  // Each processor gets NumMyElements points

  if (NumGlobal_Elements < -1) 
    throw ReportError("NumGlobal_Elements = " + toString(NumGlobal_Elements) + ".  Should be >= -1.", -1);
  if (NumMy_Elements < 0) 
    throw ReportError("NumMy_Elements = " + toString(NumMy_Elements) + ".  Should be >= 0.", -2);
  for (i = 0; i < NumMy_Elements; i++)
    if (elementSizeList[i] <= 0) 
      throw ReportError("elementSizeList["+toString(i)+"] = " + toString(elementSizeList[i]) + ". Should be > 0.", -3);
  
  BlockMapData_ = new Epetra_BlockMapData(NumGlobal_Elements, 0, indexBase, comm);
  BlockMapData_->NumMyElements_ = NumMy_Elements;
  BlockMapData_->ConstantElementSize_ = false;
  BlockMapData_->LinearMap_ = false;
  // Allocate storage for global index list and element size information

  if (NumMy_Elements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(NumMy_Elements);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
    errorcode = BlockMapData_->ElementSizeList_.Size(NumMy_Elements);
    if(errorcode != 0)
      throw ReportError("Error with ElementSizeList allocation.", -99);
  }
  // Get processor information

  int NumProc = comm.NumProc();
  
  if (NumMy_Elements > 0) {
    // Compute min/max GID and element size, number of points on this processor
    BlockMapData_->MinMyGID_ = myGlobalElements[0];
    BlockMapData_->MaxMyGID_ = myGlobalElements[0];
    BlockMapData_->MinMyElementSize_ = elementSizeList[0];
    BlockMapData_->MaxMyElementSize_ = elementSizeList[0];
    BlockMapData_->NumMyPoints_ = 0;
    for (i = 0; i < NumMy_Elements; i++) {
      BlockMapData_->MyGlobalElements_[i] = myGlobalElements[i];
      BlockMapData_->ElementSizeList_[i] = elementSizeList[i];
      BlockMapData_->MinMyGID_ = EPETRA_MIN(BlockMapData_->MinMyGID_,myGlobalElements[i]);
      BlockMapData_->MaxMyGID_ = EPETRA_MAX(BlockMapData_->MaxMyGID_,myGlobalElements[i]);
      BlockMapData_->MinMyElementSize_ = EPETRA_MIN(BlockMapData_->MinMyElementSize_,elementSizeList[i]);
      BlockMapData_->MaxMyElementSize_ = EPETRA_MAX(BlockMapData_->MaxMyElementSize_,elementSizeList[i]);
      BlockMapData_->NumMyPoints_ += elementSizeList[i];
    }
  }
  else {
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->IndexBase_ - 1;
    BlockMapData_->MinMyElementSize_ = 1;
    BlockMapData_->MaxMyElementSize_ = 1;
    BlockMapData_->NumMyPoints_ = 0;
  }

  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobal_Elements, NumMy_Elements);  

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc == 1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobal_Elements);
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumMyPoints_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->MinMyGID_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MaxMyGID_;
    BlockMapData_->MinElementSize_ = BlockMapData_->MinMyElementSize_;
    BlockMapData_->MaxElementSize_ = BlockMapData_->MaxMyElementSize_;
  }
  else if (NumProc > 1) {
    // Sum up all local element and point counts to get global counts
    int *tmp_send = new int[4];
    int *tmp_recv = new int[4];
    tmp_send[0] = BlockMapData_->NumMyElements_;
    tmp_send[1] = BlockMapData_->NumMyPoints_;
    BlockMapData_->Comm_->SumAll(tmp_send, tmp_recv, 2);
    BlockMapData_->NumGlobalElements_ =  tmp_recv[0];
    BlockMapData_->NumGlobalPoints_ = tmp_recv[1];
    
    CheckValidNGE(NumGlobal_Elements);
    
    // Use the MaxAll function to find min/max GID 
    tmp_send[0] = - BlockMapData_->MinMyGID_; // Negative signs lets us do one reduction
    tmp_send[1] =   BlockMapData_->MaxMyGID_;
    tmp_send[2] = - BlockMapData_->MinMyElementSize_;
    if (BlockMapData_->NumMyElements_ == 0) 
      tmp_send[2] = - BlockMapData_->NumGlobalPoints_; // This processor has no elements, so should not sizes.
    tmp_send[3] =   BlockMapData_->MaxMyElementSize_;
    
    BlockMapData_->Comm_->MaxAll(tmp_send, tmp_recv, 4);
    
    BlockMapData_->MinAllGID_ =      - tmp_recv[0];
    BlockMapData_->MaxAllGID_ =        tmp_recv[1];
    BlockMapData_->MinElementSize_ = - tmp_recv[2];
    BlockMapData_->MaxElementSize_ =   tmp_recv[3];
    
    delete [] tmp_send;
    delete [] tmp_recv;

    // Check for constant element size
    if (BlockMapData_->MinElementSize_==BlockMapData_->MaxElementSize_) {
      BlockMapData_->ElementSize_ = BlockMapData_->MinElementSize_;
      BlockMapData_->ConstantElementSize_ = true;
    }
    
    if (BlockMapData_->MinAllGID_ < BlockMapData_->IndexBase_)
      throw ReportError("Minimum global element index = " + toString(BlockMapData_->MinAllGID_) + 
			" is less than index base = " + toString(BlockMapData_->IndexBase_) +".", -5);
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  

  EndOfConstructorOps();
}

//==============================================================================
Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap& map)
  : Epetra_Object(map.Label()),
    BlockMapData_(map.BlockMapData_)
{
  BlockMapData_->IncrementReferenceCount();
  
  // This call appears to be unnecessary overhead.  Removed 10-Aug-2004 maherou.
  // GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.
}

//==============================================================================
bool Epetra_BlockMap::SameAs(const Epetra_BlockMap & Map) const {

  // Quickest test: See if both maps share an inner data class
  if (this->BlockMapData_ == Map.BlockMapData_) 
    return(true);


  // Next check other global properties that are easy global attributes
  if (BlockMapData_->MinAllGID_ != Map.MinAllGID() ||
      BlockMapData_->MaxAllGID_ != Map.MaxAllGID() ||
      BlockMapData_->NumGlobalElements_ != Map.NumGlobalElements() ||
      BlockMapData_->IndexBase_ != Map.IndexBase()) 
    return(false);
  
  // Last possible global check for constant element sizes
  if (BlockMapData_->ConstantElementSize_ && BlockMapData_->ElementSize_!=Map.ElementSize()) 
    return(false);

  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true
 
  int numMyElements = BlockMapData_->NumMyElements_;

  int MySameMap = 1; // Assume not needed
  
  // First check if number of element is the same in each map
  if (numMyElements != Map.NumMyElements()) MySameMap = 0;
  
  // If numMyElements is the same, check to see that list of GIDs is the same
  if (MySameMap==1) {
    if (LinearMap() && Map.LinearMap() ) {
      // For linear maps, just need to check whether lower bound is the same
      if (MinMyGID() != Map.MinMyGID() )
        MySameMap = 0;
    }
    else {
      for (int i = 0; i < numMyElements; i++) {
        if (GID(i) != Map.GID(i)) {
          MySameMap = 0;
          break;
        }
      }
    }
  }
//    for (int i = 0; i < numMyElements; i++)
//      if (GID(i) != Map.GID(i)) MySameMap = 0;

  // If GIDs are the same, check to see element sizes are the same
  if (MySameMap==1 && !BlockMapData_->ConstantElementSize_) {
    int * sizeList1 = ElementSizeList();
    int * sizeList2 = Map.ElementSizeList();
    for (int i = 0; i < numMyElements; i++) if (sizeList1[i] != sizeList2[i]) MySameMap=0;
  }
  // Now get min of MySameMap across all processors

  int GlobalSameMap = 0;
  int err = Comm().MinAll(&MySameMap, &GlobalSameMap, 1);
  assert(err==0);

  return(GlobalSameMap==1);
}

//==============================================================================
bool Epetra_BlockMap::PointSameAs(const Epetra_BlockMap & Map) const
{
  if (this->BlockMapData_ == Map.BlockMapData_) 
    return(true);
  
  if (BlockMapData_->NumGlobalPoints_ != Map.NumGlobalPoints() ) 
    return(false);
  
  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true

  int MySameMap = 1; // Assume not needed
  if (BlockMapData_->NumMyPoints_ != Map.NumMyPoints())
    MySameMap = 0;

  int GlobalSameMap = 0;
  int err = Comm().MinAll(&MySameMap, &GlobalSameMap, 1);
  assert( err == 0 );

  return(GlobalSameMap==1);
}

//==============================================================================
int Epetra_BlockMap::MyGlobalElements(int * myGlobalElements) const
{
  // If the global element list is not create, then do so.  This can only happen when
  // a linear distribution has been specified.  Thus we can easily construct the update
  // list in this case.

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;
  
  if (BlockMapData_->MyGlobalElements_.Length() == 0)
    for (i = 0; i < numMyElements; i++)
      myGlobalElements[i] = BlockMapData_->MinMyGID_ + i;
  else
    for (i = 0; i < numMyElements; i++)
      myGlobalElements[i] = BlockMapData_->MyGlobalElements_[i];
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::MyGlobalElements() const {
  int numMyElements = BlockMapData_->NumMyElements_;  

  // If ElementSizeList not built, do so
  if(BlockMapData_->MyGlobalElements_.Length() == 0 && numMyElements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(numMyElements + 1);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
    
    // Build the array
    for (int i = 0; i < numMyElements; i++)
      BlockMapData_->MyGlobalElements_[i] = BlockMapData_->MinMyGID_ + i;
  }
  return(BlockMapData_->MyGlobalElements_.Values());
}

//==============================================================================
int Epetra_BlockMap::FirstPointInElement(int lid) const
{
  if (!MyLID(lid)) 
    EPETRA_CHK_ERR(-1);
  
  int entry;

  if (ConstantElementSize())
    entry = MaxElementSize() * lid; // convert to vector entry
  else {
    int * entrylist = FirstPointInElementList(); // get entry list
    entry = entrylist[lid];
  }
  return(entry);
}

//==============================================================================
int Epetra_BlockMap::FirstPointInElementList(int * firstPointInElementList) const
{
  // If the first element entry list is not create, then do so.  

  // Note: This array is of length NumMyElement+1

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->FirstPointInElementList_.Length() == 0) {
    firstPointInElementList[0] = 0; // First element of first entry is always zero
    
    if (BlockMapData_->ConstantElementSize_)
      for (i = 0; i < numMyElements; i++)
	firstPointInElementList[i+1] = firstPointInElementList[i] + BlockMapData_->ElementSize_;
    else
      for (i = 0; i < numMyElements; i++)
	firstPointInElementList[i+1] = firstPointInElementList[i] + BlockMapData_->ElementSizeList_[i];
  }
  else 
    for (i = 0; i <= numMyElements; i++)
      firstPointInElementList[i] = BlockMapData_->FirstPointInElementList_[i];
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::FirstPointInElementList() const {
  int numMyElements = BlockMapData_->NumMyElements_;

  // If ElementSizeList not built, do so
  if ((BlockMapData_->FirstPointInElementList_.Length() == 0) && (numMyElements > 0)) {
    BlockMapData_->FirstPointInElementList_.Size(BlockMapData_->NumMyElements_ + 1);
    BlockMapData_->FirstPointInElementList_[0] = 0; // First element of first entry is always zero
    if (BlockMapData_->ConstantElementSize_)
      for (int i = 0; i < numMyElements; i++)
	BlockMapData_->FirstPointInElementList_[i+1] = BlockMapData_->FirstPointInElementList_[i] + BlockMapData_->ElementSize_;
    else
      for (int i = 0; i < numMyElements; i++)
	BlockMapData_->FirstPointInElementList_[i+1] = BlockMapData_->FirstPointInElementList_[i] + BlockMapData_->ElementSizeList_[i];
  }
  return(BlockMapData_->FirstPointInElementList_.Values());
}

//==============================================================================
int Epetra_BlockMap::ElementSizeList(int * elementSizeList) const
{
  // If the element size list is not create, then do so.  This can only happen when
  // a constant element size has been specified.  Thus we can easily construct the element size
  // list in this case.

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->ElementSizeList_.Length() == 0)
    for (i = 0; i < numMyElements; i++)
      elementSizeList[i] = BlockMapData_->ElementSize_;
  else
    for (i = 0; i < numMyElements; i++)
      elementSizeList[i] = BlockMapData_->ElementSizeList_[i];
  
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::ElementSizeList() const {
  int numMyElements = BlockMapData_->NumMyElements_;

  // If ElementSizeList not built, do so
  if ((BlockMapData_->ElementSizeList_.Length() == 0) && (numMyElements > 0)) {
    BlockMapData_->ElementSizeList_.Size(numMyElements);
    for (int i = 0; i < numMyElements; i++)
      BlockMapData_->ElementSizeList_[i] = BlockMapData_->ElementSize_;
  }
  return(BlockMapData_->ElementSizeList_.Values());
}

//==============================================================================
int Epetra_BlockMap::PointToElementList(int * pointToElementList) const {
  // Build an array such that the local element ID is stored for each point

  int i;
  if (BlockMapData_->PointToElementList_.Length() == 0) {
    int numMyElements = BlockMapData_->NumMyElements_;
    int * ptr = pointToElementList;
    for (i = 0; i < numMyElements; i++) {
      int Size = ElementSize(i);
      for (int j = 0; j < Size; j++) 
	*ptr++ = i;
    }
  }
  else {
    int numMyPoints = BlockMapData_->NumMyPoints_;
    for (i = 0; i < numMyPoints; i++)
      pointToElementList[i] = BlockMapData_->PointToElementList_[i];
  }
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::PointToElementList() const {

  // If PointToElementList not built, do so
  if ((BlockMapData_->PointToElementList_.Length() == 0) && (BlockMapData_->NumMyPoints_ > 0)) {
    BlockMapData_->PointToElementList_.Size(BlockMapData_->NumMyPoints_);
    int numMyElements = BlockMapData_->NumMyElements_;
    int * ptr = BlockMapData_->PointToElementList_.Values();
    for (int i = 0; i < numMyElements; i++) {
      int Size = ElementSize(i);
      for (int j = 0; j < Size; j++) 
	*ptr++ = i;
    }
  }
  return(BlockMapData_->PointToElementList_.Values());
}

//==============================================================================
int Epetra_BlockMap::ElementSize(int lid) const {

  if (ConstantElementSize()) 
    return(BlockMapData_->ElementSize_);
  else
    return(BlockMapData_->ElementSizeList_[lid]);
}

//==============================================================================
bool Epetra_BlockMap::IsOneToOne() const {
  if(!BlockMapData_->OneToOneIsDetermined_){
    BlockMapData_->OneToOne_ = DetermineIsOneToOne();
    BlockMapData_->OneToOneIsDetermined_ = true;
  }
  return(BlockMapData_->OneToOne_);
}

//==============================================================================
void Epetra_BlockMap::GlobalToLocalSetup()
{
  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->NumGlobalElements_ == 0) {
    return; // Nothing to do
  }

  if (LinearMap() || numMyElements == 0) {
    return; // Nothing else to do
  }

  // Build LID_ vector to make look up of local index values fast

#ifdef EPETRA_BLOCKMAP_NEW_LID

  // Here follows an optimization that checks for an initial block of
  // contiguous GIDs, and stores the GID->LID mapping for those in a
  // more efficient way.  This supports a common use case for
  // overlapping Maps, in which owned entries of a vector are ordered
  // before halo entries.
  //
  // Epetra defines EPETRA_BLOCKMAP_NEW_LID by default (see the top of
  // this file).

  //check for initial contiguous block
  int val = BlockMapData_->MyGlobalElements_[0];
  for( i = 0 ; i < numMyElements; ++i ) {
    if (val != BlockMapData_->MyGlobalElements_[i]) break;
    ++val;
  }
  BlockMapData_->LastContiguousGIDLoc_ = i - 1;
  if (BlockMapData_->LastContiguousGIDLoc_ < 0) {
    BlockMapData_->LastContiguousGID_ = BlockMapData_->MyGlobalElements_[0];
  }
  else {
    BlockMapData_->LastContiguousGID_ =
      BlockMapData_->MyGlobalElements_[BlockMapData_->LastContiguousGIDLoc_];
  }

  //Hash everything else
  if(i < numMyElements) {
    if (BlockMapData_->LIDHash_ != NULL) {
      delete BlockMapData_->LIDHash_;
    }

    BlockMapData_->LIDHash_ = new Epetra_HashTable(numMyElements - i + 1 );
    for(; i < numMyElements; ++i )
      BlockMapData_->LIDHash_->Add( BlockMapData_->MyGlobalElements_[i], i );
  }
    
#else
    
  int SpanGID = BlockMapData_->MaxMyGID_ - BlockMapData_->MinMyGID_ + 1;
  BlockMapData_->LID_.Size(SpanGID);
    
  for (i = 0; i < SpanGID; i++) 
    BlockMapData_->LID_[i] = -1; // Fill all locations with -1
    
  for (i = 0; i < numMyElements; i++) {
    int tmp = BlockMapData_->MyGlobalElements_[i] - BlockMapData_->MinMyGID_;
    assert(tmp >= 0); 
    assert(tmp < SpanGID);
    BlockMapData_->LID_[BlockMapData_->MyGlobalElements_[i] - BlockMapData_->MinMyGID_] = i; // Spread local indices
  }

#endif

}

//==============================================================================
int Epetra_BlockMap::LID(int gid) const
{
  if ((gid < BlockMapData_->MinMyGID_) || 
      (gid > BlockMapData_->MaxMyGID_)) {
    return(-1); // Out of range
  }

  if (BlockMapData_->LinearMap_) {
    return(gid - BlockMapData_->MinMyGID_); // Can compute with an offset
  }

  if( gid >= BlockMapData_->MyGlobalElements_[0] &&
      gid <= BlockMapData_->LastContiguousGID_ ) {
    return( gid - BlockMapData_->MyGlobalElements_[0] );
  }

#ifdef EPETRA_BLOCKMAP_NEW_LID
  return BlockMapData_->LIDHash_->Get( gid );
#else
  return(BlockMapData_->LID_[gid - BlockMapData_->MinMyGID_]); // Find it in LID array  
#endif
}

//==============================================================================
int Epetra_BlockMap::GID(int lid) const
{
  if ((BlockMapData_->NumMyElements_==0) ||
      (lid < BlockMapData_->MinLID_) || 
      (lid > BlockMapData_->MaxLID_)) {
    return(BlockMapData_->IndexBase_ - 1); // Out of range
  }

  if (LinearMap()) {
    return(lid + BlockMapData_->MinMyGID_); // Can compute with an offset
  }

  return(BlockMapData_->MyGlobalElements_[lid]); // Find it in MyGlobalElements array
}

//==============================================================================
int Epetra_BlockMap::FindLocalElementID(int PointID, int & ElementID, int & ElementOffset) const {

  if (PointID >= BlockMapData_->NumMyPoints_)
    return(-1); // Point is out of range
  
  if (ConstantElementSize()) {
    ElementID = PointID / BlockMapData_->MaxElementSize_;
    ElementOffset = PointID % BlockMapData_->MaxElementSize_;
    return(0);
  }
  else {
    int * tmpPointToElementList = PointToElementList();
    int * tmpFirstPointInElementList = FirstPointInElementList();
    ElementID = tmpPointToElementList[PointID];
    ElementOffset = PointID - tmpFirstPointInElementList[ElementID];
    return(0);
  }
}

//==============================================================================
int Epetra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList,
				  int * PIDList, int * LIDList,
				  int * SizeList) const
{
  if (BlockMapData_->Directory_ == NULL) {
    BlockMapData_->Directory_ = Comm().CreateDirectory(*this);
  }

  Epetra_Directory* directory = BlockMapData_->Directory_;
  if (directory == NULL) {
    return(-1);
  }

  EPETRA_CHK_ERR( directory->GetDirectoryEntries(*this, NumIDs, GIDList,
						 PIDList, LIDList, SizeList) );

  return(0);
}

//==============================================================================
bool Epetra_BlockMap::DetermineIsOneToOne() const
{
  if (Comm().NumProc() < 2) {
    return(true);
  }
  
  if (BlockMapData_->Directory_ == NULL) {
    BlockMapData_->Directory_ = Comm().CreateDirectory(*this);
  }

  Epetra_Directory* directory = BlockMapData_->Directory_;
  if (directory == NULL) {
    throw ReportError("Epetra_BlockMap::IsOneToOne ERROR, CreateDirectory failed.",-1);
  }

  return(directory->GIDsAllUniquelyOwned());
}

//==============================================================================
bool Epetra_BlockMap::IsDistributedGlobal(int numGlobalElements, int numMyElements) const {

  bool isDistributedGlobal = false; // Assume map is not global distributed
  if (BlockMapData_->Comm_->NumProc() > 1) {
    int LocalReplicated = 0;
    int AllLocalReplicated;
    if (numGlobalElements == numMyElements) 
      LocalReplicated=1;
    BlockMapData_->Comm_->MinAll(&LocalReplicated, &AllLocalReplicated, 1);
    
    // If any PE has LocalReplicated=0, then map is distributed global
    if (AllLocalReplicated != 1) 
      isDistributedGlobal = true;
  }
  return(isDistributedGlobal);
}

//==============================================================================
void Epetra_BlockMap::CheckValidNGE(int numGlobalElements) {
  // Check to see if user's value for numGlobalElements is either -1 
  // (in which case we use our computed value) or matches ours.
  if ((numGlobalElements != -1) && (numGlobalElements != BlockMapData_->NumGlobalElements_)) {
    int BmdNumGlobalElements = BlockMapData_->NumGlobalElements_;
    CleanupData();
    throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(numGlobalElements) + 
		      ".  Should equal " + toString(BmdNumGlobalElements) + 
		      ", or be set to -1 to compute automatically", -4);
  }
}

//==============================================================================
void Epetra_BlockMap::EndOfConstructorOps() {
  BlockMapData_->MinLID_ = 0;
  BlockMapData_->MaxLID_ = EPETRA_MAX(BlockMapData_->NumMyElements_ - 1, 0);
  
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.
}

//==============================================================================
void Epetra_BlockMap::Print(ostream & os) const
{
  int * MyGlobalElements1 = MyGlobalElements();
  int * FirstPointInElementList1 = 0;
  int * ElementSizeList1 = 0;
  if (!ConstantElementSize()) {
    FirstPointInElementList1 = FirstPointInElementList();
    ElementSizeList1 = ElementSizeList();
  }
  int MyPID = Comm().MyPID();
  int NumProc = Comm().NumProc();
  
  for (int iproc = 0; iproc < NumProc; iproc++) {
    if (MyPID == iproc) {
      if (MyPID == 0) {
	os <<  "\nNumber of Global Elements  = "; os << NumGlobalElements(); os << endl;
	os <<    "Number of Global Points    = "; os << NumGlobalPoints(); os << endl;
	os <<    "Maximum of all GIDs        = "; os << MaxAllGID(); os << endl;
	os <<    "Minimum of all GIDs        = "; os << MinAllGID(); os << endl;
	os <<    "Index Base                 = "; os << IndexBase(); os << endl;
	if (ConstantElementSize())
	  os <<  "Constant Element Size      = "; os << ElementSize(); os << endl;
      }
      os << endl;
      
      os <<    "Number of Local Elements   = "; os << NumMyElements(); os << endl;
      os <<    "Number of Local Points     = "; os << NumMyPoints(); os << endl;
      os <<    "Maximum of my GIDs         = "; os << MaxMyGID(); os << endl;
      os <<    "Minimum of my GIDs         = "; os << MinMyGID(); os << endl;
      os << endl;
			
      os.width(14);
      os <<  "     MyPID"; os << "    ";
      os.width(14);
      os <<  "       Local Index "; os << " ";
      os.width(14);
      os <<  "      Global Index "; os << " ";
      if (!ConstantElementSize()) {
	os.width(14);
	os <<" FirstPointInElement "; os << " ";
	os.width(14);
	os <<"   ElementSize "; os << " ";
      }
      os << endl;
			
      for (int i = 0; i < NumMyElements(); i++) {
	os.width(14);
	os <<  MyPID; os << "    ";
	os.width(14);
	os <<  i; os << "    ";
	os.width(14);
	os <<  MyGlobalElements1[i]; os << "    ";
	if (!ConstantElementSize()) {	  
	  os.width(14);
	  os << FirstPointInElementList1[i]; os << "    ";
	  os.width(14);
	  os << ElementSizeList1[i]; os << "    ";
	}
	os << endl;
      }
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
  return;
}

//==============================================================================
Epetra_BlockMap::~Epetra_BlockMap()  {
  CleanupData();
}

//==============================================================================
void Epetra_BlockMap::CleanupData()
{
  if(BlockMapData_ != 0) {

    BlockMapData_->DecrementReferenceCount();
    if(BlockMapData_->ReferenceCount() == 0) {
      delete BlockMapData_;
      BlockMapData_ = 0;
    }
  }
}

//=============================================================================
Epetra_BlockMap & Epetra_BlockMap::operator= (const Epetra_BlockMap & map)
{
  if((this != &map) && (BlockMapData_ != map.BlockMapData_)) {
    CleanupData();
    BlockMapData_ = map.BlockMapData_;
    BlockMapData_->IncrementReferenceCount();
  }

  return(*this);
}
