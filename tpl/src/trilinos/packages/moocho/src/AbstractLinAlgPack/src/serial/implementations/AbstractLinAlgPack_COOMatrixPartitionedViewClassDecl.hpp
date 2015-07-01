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

#ifndef COO_MATRIX_PARTITIONED_VIEW_CLASS_DECL_H
#define COO_MATRIX_PARTITIONED_VIEW_CLASS_DECL_H

#include "AbstractLinAlgPack_SparseCOOPtrElement.hpp"
#include "AbstractLinAlgPack_TransSparseCOOElementViewIter.hpp"
#include "MiRefCount.h"

namespace AbstractLinAlgPack {

namespace COOMatrixPartitionedViewUtilityPack {
template <class T_Indice, class T_Value> class Partition;
template <class T_Indice, class T_Value> class TransposedPartition;
}

// ///////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////
/** \brief COO Matrix partitioning view class.
  *
  * This class is used to provide a set of views (partitions)
  * of a permuted sparse coordinate (COO) matrix sorted in Fortran
  * like compatable arrays (val, ivect, jvect).
  *
  * The function #create_view(...)# or a constructor is called to setup a view.
  * In these initizliation functions the Fortran compatable COO matrix is specified as:\\
  *		(#rows#,#cols#,#nz#,#val#(#nz#),#ivect#(#nz#),#jvect#(#nz#))\\
  * The row and column permutations are given as a set of inverse permutation vectors:\\
  *		(#inv_row_perm#(#rows#),#inv_col_perm#(#cols#))\\
  * And the row and column partitioning is given as:
  *		(#num_row_part#,#row_part#(#num_row_part#),#num_col_part#,#col_part#(#num_col_part#))
  *
  * After these function execute successfully the object provides views to
  * the given partitions (singular or consecutive) of the given permuted COO matrix.
  * The partitions are ordered by column or by row depending on the value of
  * #partition_order#.
  *
  * For example, concider the following permuted and partitioned sparse COO matrix:
  \begin{verbatim}
  
         Non-permuted COO matrix		
                        
     1       2       3       4       5	
    ---     ---     ---     ---     ---	
                
1 |         1.2             1.4         |
                
2 | 2.1     2.2             2.4     2.5 |
                          
3 |                 3.3                 |
                        
4 |                                 4.5 |
                        
5 | 5.1     5.2     5.3     5.4         |
                        
6 |                 6.3             6.5 |
                        
7 | 7.1                     7.4         |
                        
8 |         8.2     8.3                 |


 ==>

      Permuted and partitioned COO matrix
 
         1       3       5       2       4
   ------------------------------------------
  |      1       2       3       4       5	
  |     ---     ---     ---     ---     ---	
  |							 				
2 | 1 | 2.1             2.5  +  2.2     2.4 |
  |                          +				
3 | 2 |         3.3          +              |
  |                          +
6 | 3 |         6.3     6.5  +              |
  |    +++++++++++++++++++++++++++++++++++++
1 | 4 |                      +  1.2     1.4 |
  |                          +
4 | 5 |                 4.5  +              |
  |                          +
5 | 6 | 5.1     5.3          +  5.2     5.4 |
  |                          +
7 | 7 | 7.1                  +          7.4 |
  |                          +
8 | 8 |         8.3          +  8.2         |
  
  \end{verbatim}
  *	The following quantities are used to specify the example shown above:
  \begin{verbatim}

rows  = 8, cols = 5, nz = 18
val   = { 2.1,5.1,7.1,1.2,2.2,5.2,8.2,3.3,5.3,6.3,8.3,1.4,2.4,5.4,7.4,2.5,4.5,6.5 }
ivect = { 2,  5,  7,  1,  2,  5,  8,  3,  5,  6,  8,  1,  2,  5,  7,  2,  4,  6   }
jvect = { 1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5   }
  
inv_row_perm  = { 4,  1,  2,  5,  6,  3,  7,  8  }
inv_col_perm  = { 1,  4,  2,  5,  3  }
  
num_row_part  = 2
row_part      = { 1,  4,  9 }
num_col_part  = 2
col_part      = { 1,  4,  6 }
  
  \end{verbatim}
  * The partitions are ordered by row or column counting from the upper left hand partition.
  * For this example the partitions are numbered as:
  \begin{verbatim}
  partition_order = PARTITION_BY_ROW          partition_order = PARTITION_BY_COL 
            --------------                             ---------------
           |   1   |  2   |                           |   1   |   3   |
           |-------|------|                           |-------|-------|
           |   3   |  4   |                           |   2   |   4   |
            --------------                             --------------- 
  \end{verbatim}
  *	The overall partition number is given by a call to #overall_part_num(row_p,col_p)# where
  * #row_p# and #col_p# are the row and column parition numbers respectively.  In the example,
  * the lower left partition has #row_p# = 2 and #col_p# = 1 and #overall_part_num(row_p,col_p)# would
  * return 3 for #partition_order = PARTITION_BY_ROW# and 2 for #partition_order = PARTITION_BY_COL#.
  *
  * Also, given the overall partition number you can obtain the row and column partition numbers
  * by calling #row_part_num(overall_p)# and #col_part_num(overall_p)# respectively.
  *
  * The row and column indices for the partitions can be extracted with calls to 
  * #get_row_part(row_part)# and #get_col_part(col_part)# respectively.
  *
  * Actual paritions or consecutive partitions are accessed by way of objects of the types
  * #partition_type# and #transposed_partition_type#.  These types conform to the 
  * COOMatrixTemplateInterface specification so they can be used with the associated templated
  * linear algebra operations.  A #partition_type# object
  * is created by a call to #partition(overall_p)# or #partition(row_p,col_p)# for
  * single partitions and #partition(Range1D(overall_p1,overall_p2))# for consecutive partitions.
  * For example you could create a partition object for the first partition by calling
  * #partition(1)# or #partition(1,1)# or #partition(Range1D(1,1))#.  You can
  * create a transposed parition object by calling the nonmember function
  * #transposed_partition_type trans(parition_type&)#.  In the above example, using
  * parition ordering by column you could create a transposed view of the two left hand
  * side partitions by calling #trans(this->partition(Range1D(1,2)))#.
  *
  * All of this gives the client a lot of flexibility in how a COO matrix is accessed.
  *
  * The default copy constructor is allowed.  The default assignment operator is not allowed.
  */
template <class T_Indice, class T_Value>
class COOMatrixPartitionedView {
public:
  
  // //////////////////////////////////////////////////////////////////////////
  /** @name Public Types */
  //@{

  /** \brief . */
  typedef COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>
                              partition_type;
  /** \brief . */
  typedef COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>
                              transposed_partition_type;
  /** \brief . */
  typedef T_Indice										indice_type;
  /** \brief . */
  typedef AbstractLinAlgPack::size_type						size_type;
  /** \brief . */
  typedef ptrdiff_t										difference_type;
  /** \brief . */
  enum EPartitionOrder { PARTITION_BY_ROW, PARTITION_BY_COL };
  /** \brief . */
  class UninitializedException: public std::logic_error
  {public: UninitializedException(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  // ///////////////////////////////////////////////////////////////////////////////////////
  /** @name Public interface */
  //@{

  /** @name Constructors and initializes
    *
    * The default copy constructor is allowed.  The constructed object will share
    * the storage space for the COO matrix view.
    */
  //@{

  /** \brief Construct with no view set.
    *
    * Postconditions:\begin{itemize}
    *	\item	#rows() == 0#
    *	\item	#cols() == 0#
    *	\item	#nz() == 0#
    *	\end{itemize}
    *
    */
  COOMatrixPartitionedView();

  /** \brief Construct with a view to a partitioned COO matrix set.
    *
    * Equivalent to calling the default constructor then #create_view(...)#.
    */
  COOMatrixPartitionedView(
        size_type				rows
      , size_type				cols
      , size_type				nz		
      , value_type			val[]
      , const indice_type		ivect[]
      , const indice_type		jvect[]
      , const size_type		inv_row_perm[]
      , const size_type		inv_col_perm[]
      , const size_type		num_row_part
      , const size_type		row_part[]
      , const size_type		num_col_part
      , const size_type		col_part[]
      , const EPartitionOrder	partition_order )

    : num_row_part_(0), num_col_part_(0)
  {
    create_view(rows, cols, nz, val, ivect, jvect, inv_row_perm, inv_col_perm, num_row_part
      , row_part, num_col_part, col_part, partition_order);
  }
    // I tried defining this function else where as both an inlined and
    // noninlined function but MS VC++ 5.0 generated
    // a linking error saying that it could not find this function.
    // ToDo: Investagate this problem further at a latter time.

  /** \brief Crete a view to a COO matrix.
    *
    * The arguments to this function specify COO matrix
    *		(#rows#,#cols#,#nz#,#val#(#nz#),#ivect#(#nz#),#jvect#(#nz#))
    * a set of row and column permutations
    *		(#inv_row_perm#(#rows#),#inv_col_perm#(#cols#))
    * and a set of row and column partitions:
    *		(#num_row_part#,#row_part#(#num_row_part#),#num_col_part#,#col_part#(#num_col_part#))
    *
    * After this function executes successfully #this# object provides views to
    * the given partitions (singular or consecutive) of the given permuted COO matrix.
    * The partitions are ordered by column or by row depending on the value of
    * #partition_order#.
    *
    * Warning!  Duplicate elements (#ivect[k1] == ivect[k2] && jvect[k1] == jvect[k2]#)
    * are not looked for in the input COO matrix.  It is up to the client to make sure
    * that there are no duplicates since this is an expensive test since the matrix
    * elements are not assumed to be sorted.
    *
    * It is important to note that #this# object keeps a pointer to the array
    * of nonzero elements in #val# passed into this function.  Therefore,
    * the numerical values of these elements can be altered and #this# object
    * will still provide the view into this coordinate matrix.  If the memory
    * pointed to by #val# change then this operation must be called again
    * setupd the view again.
    *
    * Preconditions:\begin{itemize}
    * \item #nz <= rows * cols#, (throw #std::invalid_argument#)
    * \item #rows >= ivect[#k#] >= 1#, k = 0,...,#nz-1# (throw #std::out_of_range#)
    * \item #cols >= jvect[#k#] >= 1#, k = 0,...,#nz-1# (throw #std::out_of_range#)
    * \item #num_row_part <= rows#, (throw #std::invalid_argument#)
    * \item #num_col_part <= cols#, (throw #std::invalid_argument#)
    * \item #rows >= inv_row_perm[#i#] >= 1#, i = 0,...,#rows-1# (throw #std::out_of_range#)
    * \item #cols >= inv_col_perm[#j#] >= 1#, j = 0,...,#cols-1# (throw #std::out_of_range#)
    * \item #inv_row_perm[#i1#] != inv_row_perm[#i2]#, for i1 anb i2 in set {0,...,#rows-1#}
    * \item #inv_col_perm[#j1#] != inv_col_perm[#j2]#, for j1 anb j2 in set {0,...,#cols-1#}
    * \item #row_part[#k+1#] > row_part[#k#]#, k = 0,...,#num_row_part-2# (throw #std::domain_error#)
    * \item #col_part[#k+1#] > col_part[#k#]#, k = 0,...,#num_col_part-2# (throw #std::domain_error#)
    * \item #row_part[num_row_part-1] <= rows + 1# (throw #std::out_of_range#)
    * \item #col_part[num_col_part-1] <= cols + 1# (throw #std::out_of_range#)
    * \item #partition_order == PARTITION_BY_ROW || partition_order == PARTITION_BY_COL#
    *			(throw #std::invalid_argument#)
    * \end{itemize}
    *
    * Postconditions:\begin{itemize}
    * \item ?
    * \end{itemize}
    *
    *	@param	rows	number of rows in the COO matrix.
    *	@param	cols	number of columns in the COO matrix.
    *	@param	nz		number of nonzero elements in the COO matrix.
    *	@param	val		array (length #nz#) of the values of the nonzero elements in the COO matrix.
    *	@param	ivect	array (length #nz#) of the row indices (1-based) in the original COO matrix.
    *	@param	jvect	array (length #nz#) of the column indices (1-based) in the original COO matrix.
    * @param	inv_row_perm	array (length #rows#) for the inverse row permutations.
    *							row_i_new == #inv_row_perm[#row_i_old#]#.
    * @param	inv_col_perm	array (length #cols#) for the inverse column permutations.
    *							col_j_new == #inv_col_perm[#col_j_old#]#..
    *	@param	num_row_part	number of row partitions for the partitioned view to create.
    *	@param	row_part	array (length #num_row_part#) that specifies the row partitions.
    *						The first #num_row_part# elements of #row_part# give the row indices
    *						for the starts of the partitions.  The last element gives the row
    *						indice one past the last row in the last partition.
    *	@param	num_col_part	number of column partitions for the partitioned view to create.
    *	@param	col_part	array (length #num_col_part#) that specifies the column partitions.
    *						The first #num_col_part# elements of #col_part# give the column indices
    *						for the starts of the partitions.  The last element gives the column
    *						indice one past the last column in the last partition.
    *	@param	parition_order	PARTITION_BY_ROW: the overall partitions are numbered from left to
    *							right by row starting at the upper left partition.  PARTITION_BY_COL
    *							the overall partitions are numbered from top to bottom by column
    *							starting at the upper left partition.
    */
  void create_view(
        size_type				rows
      , size_type				cols
      , size_type				nz		
      , value_type			val[]
      , const indice_type		ivect[]
      , const indice_type		jvect[]
      , const size_type		inv_row_perm[]
      , const size_type		inv_col_perm[]
      , const size_type		num_row_part
      , const size_type		row_part[]
      , const size_type		num_col_part
      , const size_type		col_part[]
      , const EPartitionOrder	partition_order );

  /** \brief Bind the view of another partitioned matrix.
    *
    * After being called #this# and #coom_view# will share storage
    * for the same COO matrix view.
    */
  void bind(const COOMatrixPartitionedView& coom_view);

  /// Free the allocated memory and make uninitialized
  void free();

  /// Returns true if a view has been initialized.
  bool is_initialized() const;

  //@}

  /** @name Partitioning information
    *
    * All of these member functions have the followin preconditions:
    *
    * Preconditions:\begin{itemize}
    *	\item #is_initialized() == true# (throw #UninitializedException#)
    * \end{itemize}
    * */
  //@{

  /// return the number of rows of the total view
  size_type rows() const;

  /// return the number of columns in the total view
  size_type cols() const;

  /// return the number of nonzero elements in the total view
  size_type nz() const;

  /// return the number row partitions
  size_type num_row_part() const;
  
  /// return the number column partitions
  size_type num_col_part() const;

  /** \brief get the array of row partitions.
    *
    *	@param	row_part	[O]	array (length #this->num_row_part() + 1#) of
    *							row indices of start of the row partitions.
    */
  void get_row_part(indice_type row_part[]) const;

  /** \brief get the array of column partitions.
    *
    *	@param	col_part	[O]	array (length #this->num_col_part() + 1#) of
    *							column indices of start of the column partitions.
    */
  void get_col_part(indice_type col_part[]) const;

  /// Returns whether the paritions or sorted by row or by column.
  EPartitionOrder partition_order() const;

  /** \brief Returns the overall partition number (1 to #num_row_part() * num_col_part()#)
    * given the row (1 to #num_row_part()#) and column (1 to #num_col_part()#) partition numbers.
    */
  size_type overall_part_num(size_type row_p, size_type col_p) const;

  /// Returns the row parition number (1-based) given the overall partition number (1-based)
  size_type row_part_num(size_type overall_p) const;

  /// Returns the column parition number (1-based) given the overall partition number (1-based)
  size_type col_part_num(size_type overall_p) const;

  //@}

  /** @name Partition access
    *
    * All of these member functions have the followin preconditions:
    *
    * Preconditions:\begin{itemize}
    *	\item #is_initialized() == true# (throw #UninitializedException#)
    * \end{itemize}
    */
  //@{

//	///
//	/** Allow an implicit conversion from a COOMatrixPartitionedView
//	  * to a partition_type object.
//	  *
//	  * This conversion is equivalent to calling #partition(Range1D())#.
//	  */
//	operator partition_type();
//
//	///
//	operator const partition_type() const;

  /** \brief Return a partition object for the entire view.
    */
  partition_type operator()();

  /** \brief . */
  const partition_type operator()() const;

  /** \brief Return a partition object for a parition given its overall partition number (1-based).
    */
  partition_type partition(size_type overall_p);

  /** \brief . */
  const partition_type partition(size_type overall_p) const;

  /** \brief Return a partition object for a parition given its row and column
    * partition numbers.
    */
  partition_type partition(size_type row_p, size_type col_p);

  /** \brief . */
  const partition_type partition(size_type row_p, size_type col_p) const;

  /** \brief Return a partition object for a range of continous overall partition numbers.
    *
    * This is probably only usefull when a set of rows of submatrices or column
    * submatrices need to be accessed.  For example, if the paritions are
    * ordered by column and if you want to access the first n columns of partition
    * submatrices you would call.
    *
    * #this->partition(Range1D(1,this->overall_part_num(this->num_row_part(),n)))#
    *
    *	@param	rng_overall_p	[I]	the range of overall partitions.
    *								If #rng_overall_p.full_range() == true# then
    *								all of the partitions will be returned.
    */
  partition_type partition(Range1D rng_overall_p);

  /** \brief . */
  const partition_type partition(Range1D rng_overall_p) const;

  //@}

  //@}

private:
  // ////////////////////////////////////////////////////////////////
  // Private types

  typedef SparseCOOPtrElement<T_Indice,T_Value>	element_type;
  typedef std::vector<indice_type>				vector_indice_type;
  typedef std::vector<size_type>					vector_size_type;
  typedef std::vector<element_type>				ele_type;
  typedef MemMngPack::RefCount<
        vector_indice_type>					ref_vector_indice_type;
  typedef MemMngPack::RefCount<
        vector_size_type>					ref_vector_size_type;
  typedef MemMngPack::RefCount<
        ele_type>							ref_ele_type;

  // ///////////////////////////////////////////////////////////////
  // Private data members

  size_type					num_row_part_;
  // The number of partions the COO matrix is sliced up into by rows.

  size_type					num_col_part_;
  // The number of partions the COO matrix is sliced up into by columns.

  ref_vector_size_type		ref_row_part_;
  //
  // row_part = { row_p_1_i, row_p_2_i,...,row_p_nrp_i, row_p_nrp_last_i + 1 }
  //		where:	row_p_k_i is the row indice for row partition k
  //				nrp = num_row_part
  //				row_p_nrp_last_i is the last row indice in the last row partition
  //
  // vector (length num_row_part_ + 1) of row indices for the start of the row partitions.
  // row_part[i - 1] is the row indice for the start of the ith row partition (i = 1,...,num_row_part_).
  // The last indice gives the row indice that is one past the last row of the last
  // row partition.  It is included for the efficient calculation of the number of rows
  // in a partition or a set of partitions.

  ref_vector_size_type		ref_col_part_;
  //
  // col_part = { col_p_1_i, col_p_2_i,...,col_p_nrp_i, col_p_nrp_last_i + 1 }
  //		where:	col_p_k_i is the col indice for col partition k
  //				nrp = num_col_part
  //				col_p_nrp_last_i is the last col indice in the last col partition
  //
  // vector (length num_row_part_ + 1) of column indices for the start of the columns partitions.
  // col_part[j - 1] is the column indice for the start of the jth column partition (j = 1,...,num_col_part_).
  // The last indice gives the column indice that is one past the last column of the last
  // column partition.  It is included for the efficient calculation of the number of columns
  // in a partition or a set of partitions.
  
  EPartitionOrder				partition_order_;
  // Specifies wheather the partitions are to be orded by columns or rows

  ref_ele_type				ref_ele_;
  // ToDo: replace with a RefCount<std::vector<element_type>> when you get compiler bug fix.
  // The storage array (dynamically allocated) for the elements in the COO matrix.
  // The elements in this array are ordered into partitions defined by the partition_order_.
  //
  //				ele_
  //		---------------------
  //		|pval	row_i	col_i| 
  //		||  |	|   |	|   ||-
  //		||  |	|   |	|   ||	parition	1
  //		||  |	|   |	|   ||
  //		||  |	|   |	|   ||-
  //		||  |	|   |	|   ||	parition	2
  //		||  |	|   |	|   ||
  //		||  |	|   |	|   ||-
  //		||  |	|   |	|   ||	.
  //		||  |	|   |	|   ||	.
  //		||  |	|   |	|   ||	.
  //		||  |	|   |	|   ||
  //		||  |	|   |	|   ||-
  //		||  |	|   |	|   ||	partition	num_row_part * num_col_part
  //		||  |	|   |	|   ||
  //		||  |	|   |	|   ||-
  //		---------------------

  ref_vector_size_type		ref_part_start_;
  //
  // part_start = { start_1 = 0, start_2,..., start_n, total_nz }
  //
  // vector (length num_row_part * num_col_part_ + 1) of the start for the elements in
  // each partition.  ele_[part_start_[overall_p - 1]] gives the first element of partion
  // overal_p and ele_[part_start_[overall_p] - 1] gives the last element in the partition.
  // This array is also used to calculate the number of nonzero elements in each parition
  // and in the overall partition.  The last element part_start_[num_row_part * num_col_part_]
  // gives the total number of nonzeros in all of the partitions and is also used
  // to allow the calculation of the number of nonzero elements in the last partition.

  // When the view is uninitialized then ele_, row_part_, col_part_ and part_start_ will
  // all uninitialized.

  // ///////////////////////////////////////////////////////////////
  // Private member functions

  // assert that the object is initialized
  void assert_initialized() const;

  size_type imp_overall_part_num(size_type row_p, size_type col_p) const;

  size_type imp_row_part_num(size_type overall_p) const;

  size_type imp_col_part_num(size_type overall_p) const;

  // Return the partition number given an indice.
  // Passing the partition vector is an optimization
  // that allows us to only call .obj() once
  // on the reference object and thus same some
  // unnecessary work.
  size_type part_num(const vector_size_type& part, size_type indice);

  // Return the overall partition number
  size_type overall_p_from_ij(const vector_size_type& row_part
    , const vector_size_type& col_part, size_type i, size_type j);

  // Return a non-const partition.  This is the function
  // that implements all partition creations.
  partition_type create_partition(Range1D rng_overall_p) const;

  // Not defined and not to be called
  COOMatrixPartitionedView& operator=(const COOMatrixPartitionedView&);

};	// end class COOMatrixPartitionedView


namespace COOMatrixPartitionedViewUtilityPack {

/** @name COOMatrixPartitionedViewUtilityPack
  * @memo C++ namespace
  */
//@{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
/** \brief Class for a partition or a set of continous
  * partitions in a partitioned COO matrix.
  *
  * This class represents the abstraction of the
  * submatrx or submatrices in a partion or
  * a set of continous partitions in a COO matrix.
  *
  * Its interface conforms to the template specification
  * COOMatrixTemplateInterface so that it can be used
  * with all of the linear algebra functions defined
  * for that interface.
  */
template <class T_Indice, class T_Value>
class Partition {
public:
  // /////////////////////////////////////////////////////
  /** @name Public types */
  //@{
  
  /** \brief . */
  typedef AbstractLinAlgPack::size_type				size_type;
  /** \brief . */
  typedef Partition<T_Indice,T_Value>				partition_type;
  /** \brief . */
  typedef ptrdiff_t								difference_type;
  /** \brief . */
  typedef SparseCOOPtrElement<T_Indice,T_Value>	element_type;
  /** \brief . */
  typedef element_type*							iterator;
  /** \brief . */
  typedef const element_type*						const_iterator;

  //@}

  // /////////////////////////////////////////////////////
  /** @name Public interface */
  //@{

  /** @name Constructors and initializes
    *
    * The default copy constructor is allowed since it has the
    * proper sematics.
    */
  //@{
  
  /** \brief Constructs an uninitialized partition view.
    *
    * Postconditions:\begin{itemize}
    * \item #rows() == 0#
    * \item #cols() == 0#
    * \item #nz() == 0#
    * \item #row_offset() == 0#
    * \item #col_offset() == 0#
    * \end{itemize}
    *
    */
  Partition();

  /** \brief Construct with the COO matrix initialized.
    *
    * Equivalent to calling the default constructor and #initialize(...)#.
    */
  Partition(
        size_type			rows
      , size_type			cols
      , size_type			nz
      , element_type*		ele
      , difference_type	row_offset
      , difference_type	col_offset );

  /** \brief Initialize the COO matrix.
    *
    * @param	rows	number of rows in the submatrix
    *	@param	cols	number of columns in the submatrix
    *	@param	nz		number of nonzero elements in the submatrix
    *	@param	ele		pointer to the array of nonzero elements in the submatrix
    *	@param	row_offset	offset for each nonzero row indice in #ele#.
    *						i = #ele#[...].row_i + #row_offset#
    *	@param	row_offset	offset for each nonzero column indice in #ele#.
    *						j = #ele#[...].col_j + #col_offset#
    */
  void initialize(
        size_type			rows
      , size_type			cols
      , size_type			nz
      , element_type*		ele
      , difference_type	row_offset
      , difference_type	col_offset );

  /** \brief bind to a partion.
    *
    * ToDo:  finish documentation for this function
    */
  void bind(const partition_type& partition);

  //@}

  /** @name COOMatrixTemplateInterface interface */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const ;
  /** \brief . */
  size_type nz() const;
  /** \brief . */
  difference_type row_offset() const;
  /** \brief . */
  difference_type col_offset() const;
  /** \brief . */
  iterator begin();
  /** \brief . */
  const_iterator begin() const;
  /** \brief . */
  iterator end();
  /** \brief . */
  const_iterator end() const;

  //@}

  //@}

private:
  // //////////////////////////////////////////////////////
  // Private types

  // //////////////////////////////////////////////////////
  // Private data members

  size_type		rows_,	// The number of rows in this COO matrix
          cols_,	// The number of columns in this COO matrix
          nz_;	// The number of nonzero elements in this COO matrix
  element_type	*ele_;	// pointer to array of elements in this COO matrix
  difference_type	row_offset_,	// offset for each row indice stored in ele
          col_offset_;	// offset for each column indice stored in ele

  // //////////////////////////////////////////////////////
  // Private member functions

  // assert that we are initialized
  void assert_initialized() const;

  // not defined and not to be called
  Partition& operator=(const Partition&);

};	// end class Partition

// /////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////
/** \brief Class for the transpose of a Partition
  *
  * This class represents the abstraction of the
  * transpose of the matrix given by a Partition<> class.
  *
  * Its interface conforms to the template specification
  * COOMatrixTemplateInterface so that it can be used
  * with all of the linear algebra functions defined
  * for this interface.
  *
  * This is a very light weight class and is just
  * as efficient as Partition when used with templated
  * linear algebra operations.
  *
  * The default constructor is allowed in which case it is initialized
  * to a default constructed partition object.  Also to the default
  * copy constructor is allowed.  The default assignment operator however
  * is not allowed since its meaning is a little confusing.  When a client
  * wants to copy the underlying partition, it should use the #bind()#
  * member function instead.
  */
template <class T_Indice, class T_Value>
class TransposedPartition {
public:
  // /////////////////////////////////////////////////////
  /** @name Public types */
  //@{
  
  /** \brief . */
  typedef Partition<T_Indice,T_Value>				partition_type;
  /** \brief . */
  typedef AbstractLinAlgPack::size_type				size_type;
  /** \brief . */
  typedef ptrdiff_t								difference_type;
  /** \brief . */
  typedef SparseCOOPtrElement<T_Indice,T_Value>	element_type;
  /** \brief . */
  typedef TransSparseCOOElementViewIter<
          element_type*
        , std::random_access_iterator_tag
        , typename element_type::indice_type
        , typename element_type::value_type&
        , difference_type>					iterator;
  /** \brief . */
  typedef TransSparseCOOElementViewIter<
          const element_type*
        , std::random_access_iterator_tag
        , typename element_type::indice_type
        , const typename element_type::value_type&
        , difference_type>					const_iterator;

  //@}

  // /////////////////////////////////////////////////////
  /** @name Public interface */
  //@{

  /** @name Constructors and initializes
    *
    * The default copy constructor is allowed since it has the
    * proper sematics.
    */
  //@{

  /** \brief Construct with the partition initialized.
    *
    * ToDo:  finish documentation for this function
    */
  TransposedPartition(const partition_type& partition);

  /** \brief bind to a partion.
    *
    * ToDo:  finish documentation for this function
    */
  void bind(const partition_type& partition);

  //@}

  /** @name COOMatrixTemplateInterface interface */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  /** \brief . */
  difference_type row_offset() const;
  /** \brief . */
  difference_type col_offset() const;
  /** \brief . */
  iterator begin();
  /** \brief . */
  const_iterator begin() const;
  /** \brief . */
  iterator end();
  /** \brief . */
  const_iterator end() const;

  //@}

  //@}

private:
  partition_type			partition_;	// actually stores a partition object.

  // Not defined and not to be called
  TransposedPartition&	operator=(const TransposedPartition&);

};	// end class TransposedPartition

//	end COOMatrixPartitionedViewUtilityPack
//@}

}	// end namespace COOMatrixPartitionedViewUtilityPack


// //////////////////////////////////////////////////////////////////////////
// Nonmember functions

/** \brief Create a transposed view of a partition object.
  */
template<class T_Indice, class T_Value>
inline COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>
trans(COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& part)
{
  typedef COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>
        transposed_partition_type;
  return transposed_partition_type(part);
}

// ///////////////////////////////////////////////////////////////////////////
// Inline member function definitions

namespace COOMatrixPartitionedViewUtilityPack {

// ///////////////////////////////////////////////////////////////////////////
// Inline members for class COOMatrixPartitionedViewUtilityPack::Partition<>

// Constructors and initializes

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::Partition()
  : rows_(0), cols_(0), nz_(0), ele_(0), row_offset_(0), col_offset_(0)
{}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::Partition(
      size_type			rows
    , size_type			cols
    , size_type			nz
    , element_type*		ele
    , difference_type	row_offset
    , difference_type	col_offset )
  : rows_(rows), cols_(cols), nz_(nz), ele_(ele), row_offset_(row_offset)
    , col_offset_(col_offset)
{}

template <class T_Indice, class T_Value>
inline void Partition<T_Indice,T_Value>::initialize(
      size_type			rows
    , size_type			cols
    , size_type			nz
    , element_type*		ele
    , difference_type	row_offset
    , difference_type	col_offset )
{
  rows_		= rows;
  cols_		= cols;
  nz_			= nz;
  ele_		= ele;
  row_offset_	= row_offset;
  col_offset_	= col_offset;
}

template <class T_Indice, class T_Value>
inline void Partition<T_Indice,T_Value>::bind(const partition_type& partition) {
  initialize(partition.rows_,partition.cols_,partition.nz_,partition.ele_
    ,partition.row_offset_,partition.col_offset_);
}

// COOMatrixTemplateInterface interface 

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::size_type
Partition<T_Indice,T_Value>::rows() const
{
  return rows_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::size_type
Partition<T_Indice,T_Value>::cols() const
{
  return cols_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::size_type
Partition<T_Indice,T_Value>::nz() const
{
  return nz_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::difference_type
Partition<T_Indice,T_Value>::row_offset() const
{
  return row_offset_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::difference_type
Partition<T_Indice,T_Value>::col_offset() const
{
  return col_offset_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::iterator
Partition<T_Indice,T_Value>::begin()
{
  assert_initialized();
  return ele_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::const_iterator
Partition<T_Indice,T_Value>::begin() const
{
  assert_initialized();
  return ele_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::iterator
Partition<T_Indice,T_Value>::end()
{
  assert_initialized();
  return ele_ + nz_;
}

template <class T_Indice, class T_Value>
inline Partition<T_Indice,T_Value>::const_iterator
Partition<T_Indice,T_Value>::end() const
{
  assert_initialized();
  return ele_ + nz_;
}

// Private member functions

template <class T_Indice, class T_Value>
inline void Partition<T_Indice,T_Value>::assert_initialized() const {
  if(!ele_)
    throw std::logic_error("Partition<...> :"
        "The COO matrix was not initizlized");
}


// ///////////////////////////////////////////////////////////////////////////
// Inline members for class COOMatrixPartitionedViewUtilityPack::TransposedPartition<>

// Constructors and initializes

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::TransposedPartition(const partition_type& partition)
  : partition_(partition)
{}

template <class T_Indice, class T_Value>
inline void TransposedPartition<T_Indice,T_Value>::bind(const partition_type& partition) {
  partition_.bind(partition);
}

// COOMatrixTemplateInterface interface

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::size_type
TransposedPartition<T_Indice,T_Value>::rows() const
{
  return partition_.cols();
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::size_type
TransposedPartition<T_Indice,T_Value>::cols() const
{
  return partition_.rows();
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::size_type
TransposedPartition<T_Indice,T_Value>::nz() const
{
  return partition_.nz();
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::difference_type
TransposedPartition<T_Indice,T_Value>::row_offset() const
{
  return partition_.col_offset();
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::difference_type
TransposedPartition<T_Indice,T_Value>::col_offset() const
{
  return partition_.row_offset();
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::iterator
TransposedPartition<T_Indice,T_Value>::begin()
{
  return iterator(partition_.begin());
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::const_iterator
TransposedPartition<T_Indice,T_Value>::begin() const
{
  return const_iterator(partition_.begin());
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::iterator
TransposedPartition<T_Indice,T_Value>::end()
{
  return iterator(partition_.end());
}

template <class T_Indice, class T_Value>
inline TransposedPartition<T_Indice,T_Value>::const_iterator
TransposedPartition<T_Indice,T_Value>::end() const
{
  return const_iterator(partition_.end());
}

}	// end namespace COOMatrixPartitionViewUtilityPack


// ///////////////////////////////////////////////////////////////////////////
// Inline members for class COOMatrixPartitionedView<>

// Constructors and initializes

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::COOMatrixPartitionedView()
  : num_row_part_(0), num_col_part_(0)
{}

template <class T_Indice, class T_Value>
inline bool COOMatrixPartitionedView<T_Indice,T_Value>::is_initialized() const {
  return ref_ele_.has_ref_set();
}

// Partitioning information

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::rows() const
{
  assert_initialized();
  const std::vector<size_type> &row_part = ref_row_part_.const_obj();
  return row_part[num_row_part_] - row_part[0];
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::cols() const
{
  assert_initialized();
  const std::vector<size_type> &col_part = ref_col_part_.const_obj();
  return col_part[num_col_part_] - col_part[0];
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::nz() const
{
  assert_initialized();
  return ref_part_start_.const_obj()[num_row_part_ * num_col_part_];
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::num_row_part() const
{
  return num_row_part_;
}
  
template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::num_col_part() const
{
  return num_col_part_;
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::EPartitionOrder
COOMatrixPartitionedView<T_Indice,T_Value>::partition_order() const
{
  assert_initialized();
  return partition_order_;
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::overall_part_num(size_type row_p, size_type col_p) const
{
  assert_initialized();
  return imp_overall_part_num(row_p, col_p);
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::row_part_num(size_type overall_p) const
{
  assert_initialized();
  return imp_row_part_num(overall_p); 
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::col_part_num(size_type overall_p) const
{
  assert_initialized();
  return imp_col_part_num(overall_p); 
}

// Partition access

//template <class T_Indice, class T_Value>
//inline COOMatrixPartitionedView<T_Indice,T_Value>::operator
//COOMatrixPartitionedView<T_Indice,T_Value>::partition_type()
//{
//	return partition(Range1D());
//}

//template <class T_Indice, class T_Value>
//inline COOMatrixPartitionedView<T_Indice,T_Value>::operator
//const COOMatrixPartitionedView<T_Indice,T_Value>::partition_type() const
//{
//	return partition(Range1D());
//}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::operator()()
{
  return partition(Range1D());
}

template <class T_Indice, class T_Value>
inline const COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::operator()() const
{
  return partition(Range1D());
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(size_type overall_p)
{
  return partition(Range1D(overall_p,overall_p));
}

template <class T_Indice, class T_Value>
inline const COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(size_type overall_p) const
{
  return partition(Range1D(overall_p,overall_p));
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(size_type row_p, size_type col_p)
{
  return partition(overall_part_num(row_p,col_p));
}

template <class T_Indice, class T_Value>
inline const COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(size_type row_p, size_type col_p) const
{
  return partition(overall_part_num(row_p,col_p));
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(Range1D rng_overall_p)
{
  return create_partition(rng_overall_p);
}

template <class T_Indice, class T_Value>
inline const COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::partition(Range1D rng_overall_p) const
{
  return create_partition(rng_overall_p);
}

// Private member functions

template <class T_Indice, class T_Value>
inline void COOMatrixPartitionedView<T_Indice,T_Value>::assert_initialized() const {
  if(!is_initialized())
    throw UninitializedException("COOMatrixPartitionedView<..>::assert_initialized() :"
        " The partitioned view has not been initialized.");
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::imp_overall_part_num(size_type row_p, size_type col_p) const
{
  return (partition_order_ == PARTITION_BY_ROW) ?
        (row_p - 1) * num_col_part_ + col_p :
        (col_p - 1) * num_row_part_ + row_p; 
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::imp_row_part_num(size_type overall_p) const
{
  return (partition_order_ == PARTITION_BY_ROW) ?
        (overall_p - 1) / num_col_part_ + 1 :
        (overall_p - 1) % num_row_part_ + 1 ; 
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::imp_col_part_num(size_type overall_p) const
{
  return (partition_order_ == PARTITION_BY_COL) ?
        (overall_p - 1) / num_row_part_ + 1 :
        (overall_p - 1) % num_col_part_ + 1 ; 
}

template <class T_Indice, class T_Value>
inline COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::overall_p_from_ij(const vector_size_type& row_part
  , const vector_size_type& col_part, size_type i, size_type j)
{
  return imp_overall_part_num(   part_num(row_part,i)
                 , part_num(col_part,j) );
}

} // end namespace AbstractLinAlgPack 

#endif // COO_MATRIX_PARTITIONED_VIEW_CLASS_DECL_H
