/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Matrix-Format dependent Functions                                    */
/* ******************************************************************** */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_mat_formats.h"
#include "ml_memory.h"
/*********************************************************************/
/* Perform a matrix vector product with an MSR matrix: c = matrix * b*/
/* NOTE: matrix can only have pre-communication and both 'recv_list' */
/*       and the 'remap' vectors within 'pre_comm' must be NULL.     */
/*                                                                   */
/* PARAMETERS                                                        */
/* ==========                                                        */
/*   matrix     On input, sparse matrix.                             */
/*   Nrows      On input, number of rows in matrix.                  */
/*   b          On input, vector to be multiply with matrix.         */
/*   c          On output, c = matrix * b.                           */
/*   Nsend      On input, number of messages sent to update vector.  */
/*********************************************************************/

void ML_restricted_MSR_mult(ML_Operator *matrix, int Nrows,
                      double b[], double c[], int Nsend)
{
   int i, *bindx, bindx_row, nzeros, j, k;
   double *val;
   struct ML_CSR_MSRdata *temp;
   ML_Comm *comm;

   if (Nrows == -57)
     ML_avoid_unused_param( (void *) &Nsend);
   comm = matrix->comm;
   ML_exchange_bdry(b, matrix->getrow->pre_comm, Nrows, comm,
                    ML_OVERWRITE,NULL);
   temp  = (struct ML_CSR_MSRdata *) matrix->data;
   val   = temp->values;
   bindx = temp->columns;

   for (i = 0; i< Nrows; i++) {

     /* compute diagonal contribution */

     *c = val[i] * b[i];

     /* nonzero off diagonal contribution */

     bindx_row = bindx[i];
     nzeros    = bindx[i+1] - bindx_row;

     for (j = 0; j < nzeros; j++) {
       k   = bindx_row + j;
       *c += val[k] * b[bindx[k]];
     }
     c++;
   }
}

/*********************************************************************/
/* Free a 'struct ML_CSR_MSRdata' associated with the matrix->data field */
/*********************************************************************/

void ML_RECUR_CSR_MSRdata_Destroy(ML_Operator *matrix)
{
   if (matrix->sub_matrix != NULL)
      ML_RECUR_CSR_MSRdata_Destroy(matrix->sub_matrix);
   ML_CSR_MSRdata_Destroy(matrix->data);
}
void ML_OnlyFreeTopLevelDataPtr(void *data)
{
  struct ML_CSR_MSRdata *temp;
  temp = (struct ML_CSR_MSRdata *) data;

  if (temp != NULL) {
    ML_free(temp);
  }
}

void ML_CSR_MSRdata_Destroy(void *data)
{
   struct ML_CSR_MSRdata *temp;

   temp = (struct ML_CSR_MSRdata *) data;
   if (temp != NULL) {
      if (temp->columns != NULL) ML_free(temp->columns);
      if (temp->values  != NULL) ML_free(temp->values);
      if (temp->rowptr  != NULL) ML_free(temp->rowptr);
      ML_free(temp);
   }
}

void ML_CSR_MSRdata_Destroy_StructOnly(void *data)
{
   struct ML_CSR_MSRdata *temp;

   temp = (struct ML_CSR_MSRdata *) data;
   if (temp != NULL) {
      ML_free(temp);
   }
}
/*********************************************************************/
/* Free a 'struct ML_VBRdata' associated with the matrix->data field */
/*********************************************************************/

void ML_RECUR_VBRdata_Destroy(void *data)
{
   ML_Operator *matrix = (ML_Operator*) data;
   if (matrix->sub_matrix != NULL)
      ML_VBRdata_Destroy(matrix->sub_matrix);
   ML_VBRdata_Destroy(matrix->data);
}

void ML_VBRdata_Destroy(void *data)
{
   struct ML_vbrdata *temp;

   temp = (struct ML_vbrdata *) data;
   if (temp != NULL) {
      if (temp->bindx != NULL) ML_free(temp->bindx);
      if (temp->val  != NULL) ML_free(temp->val);
      if (temp->rpntr  != NULL) ML_free(temp->rpntr);
      if (temp->cpntr  != NULL) ML_free(temp->cpntr);
      if (temp->bpntr  != NULL) ML_free(temp->bpntr);
      if (temp->indx  != NULL) ML_free(temp->indx);
      ML_free(temp);
   }
}
/*********************************************************************/
/* Scale the rows of the generalized CSR matrix using scale_factors. */
/*********************************************************************/

void ML_Scale_CSR(ML_Operator *input_matrix, double scale_factors[],
                  int mult_or_divide)
{
   int         i, row;
   ML_Operator *next, *current;
   struct ML_CSR_MSRdata *tmatrix;
   int    *rowptr, actual_row;
   double *values, dtemp;

   for (row = 0; row < input_matrix->getrow->Nrows; row++) {

      if (mult_or_divide == 0)
         dtemp = 1./scale_factors[row];
      else dtemp = scale_factors[row];

      actual_row = row;
      if (input_matrix->getrow->row_map != NULL)
         actual_row = input_matrix->getrow->row_map[row];


      if (actual_row != -1) {

         current = input_matrix;
         next    = current->sub_matrix;

         while ( (next != NULL) && (actual_row < next->getrow->Nrows) ) {
            current = next;
            next = next->sub_matrix;
         }
         if (next != NULL) actual_row -= next->getrow->Nrows;

         tmatrix = (struct ML_CSR_MSRdata *) current->data;
         values = tmatrix->values;
         rowptr = tmatrix->rowptr;
         for (i = rowptr[actual_row] ; i < rowptr[actual_row+1] ; i++)
            values[i] *= dtemp;
      }
   }
}

/*********************************************************************/
/* Get one block matrix row(requested_row)                           */
/* from the user's matrix and return this information  in            */
/* 'blocks, bindx, indx, values, row_length, blocks'.                */
/* If there is not enough space to complete this operation, allocate */
/* it.  Otherwise, return 1.                                         */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/* data             On input, points to user's data containing       */
/*                  matrix values.                                   */
/* requested_row    On input, requested_row                          */
/*                  gives the row indice of the block row for which  */
/*                  nonzero blocks are returned.                     */
/* row_length       On output, row_length is the number of           */
/*                  elements in the nonzero blocks in the row        */
/*                  'requested_row'                                  */
/* blocks           On output, blocks is the number of nonzero blocks*/
/*                  in the row 'requested_row'                       */
/* bindx, indx, values                                               */
/*                  On output, bindx[columns] is the number of c and */
/*                  values[k] contains the                           */
/*                  column number and value of a matrix nonzero where*/
/*                  all the nonzeros for requested_rows[0] appear    */
/*                  first followed by the nonzeros for               */
/*                  requested_rows[1], etc. NOTE: these arrays are   */
/*                  of size 'allocated_space'.                       */
/*                  Blocksize[k] contains the number of elements in  */
/*                  each non-zero block                              */
/* int_space, values_space                                           */
/*                  On input, int_space indicates the space available*/
/*                  in 'bindx', 'indx' and values_space the space in */
/*                  'values' for storing nonzeros.  If more space is */
/*                  needed, it is allocated.  The size of the arrays */
/*                  are returned on output                           */
/* index, index2    If index and index are not zero then tack the row*/
/*                  being retireved onto the end of the bindx, indx  */
/*                  and values arrays.  index is used for bindx and  */
/*                  indx while index2 is used for values.            */
/*********************************************************************/

/*Another version of this function exists without the index and with
  fewer values passed by reference in version 1.73 if this file and
  associated checkins*/

int VBR_block_getrow(ML_Operator *data, int requested_row,
   int *int_space, int *values_space, int *blocks, int **bindx,
   int **indx, double **values, int *row_length, int index, int index2)
{
   struct ML_vbrdata *input_matrix;
   int  i, bpntr_rows, bpntr_rows_1;
   int *t1, *t2, *bindx_ptr, *indx_ptr;
   int *bindx_old, *indx_old, *bpntr, start, finish;
   double *val, *t3, *values_ptr;
   int done = 0;
   int col_size;
   ML_Operator *place_holder;


#ifdef DEBUG2
   if ( (requested_row >= input_matrix->getrow->N_block_rows) || (requested_row < 0) ) {
      row_length = 0;
      return;
   }
#endif

   /*cathing a potential problem but not sure if this is a problem or would ever happen*/
   /*if (data->getrow->row_map != NULL)
   {
     printf("It is unclear as to what a non-NULL rowmap means in a blockrow matrix.\n");
     exit(1);
   }*/
   /*it is unclear what this does especially in the block case so it is commented out and trapped as an error*/
   /*if (data->getrow->row_map != NULL) {
      if (data->getrow->row_map[requested_row] != -1)
         requested_row = data->getrow->row_map[requested_row];
      else { row_length = 0;
    ML_avoid_unused_param( (void *) &row_length);
    return;}
   }*/


   place_holder = data;
   /*Find matrix to query*/
   while(!done)
   {
     if(place_holder->sub_matrix == NULL)
       done = 1;
     else
     {
       if(requested_row < place_holder->sub_matrix->getrow->N_block_rows)
         place_holder = place_holder->sub_matrix;
       else
       {
         requested_row -= place_holder->sub_matrix->getrow->N_block_rows;
         done = 1;
       }
     }
   }

   /*set alliases*/
   input_matrix = (struct ML_vbrdata *) ML_Get_MyGetrowData(place_holder);
   bpntr  = input_matrix->bpntr;
   indx_old = input_matrix->indx;

   /*store frquently used data for faster access*/
   bpntr_rows = bpntr[requested_row];
   bpntr_rows_1 = bpntr[requested_row+1];
   start = indx_old[bpntr_rows];
   finish = indx_old[bpntr_rows_1];

   /*calculate block row information*/
   *row_length = finish - start;
   *blocks = bpntr_rows_1 - bpntr_rows;

   /*Make sure we have enough space to store the block row variables we need
     one more space for the last spot of the indx*/
   if(*blocks + index >= *int_space){
     *int_space = 2*(*int_space) +1;
     if(*blocks + index >= *int_space)\
        *int_space = *int_space + 5 + index;
      t1 = (int    *) ML_allocate(*int_space*sizeof(int   ));
      t2 = (int    *) ML_allocate(*int_space*sizeof(int   ));
      if (t2 == NULL) {
         printf("Not enough space to get a matrix blockrow. A row length of \n");
         printf("%d was not sufficient\n",(*int_space-1)/2);
         fflush(stdout);
         exit(1);
      }
      for (i = 0; i < index; i++) t1[i] = (*indx)[i];
      for (i = 0; i < index; i++) t2[i] = (*bindx)[i];
      ML_free(*indx); ML_free(*bindx);
      *indx = t1;
      *bindx = t2;
   }
   if (*row_length + index2 > *values_space) {
     *values_space = 2*(*values_space) +1;
     if(*row_length + index2 > *values_space)\
        *values_space = *values_space + 5 + index2;
      t3 = (double *) ML_allocate(*values_space*sizeof(double));
      if (t3 == NULL) {
         printf("Not enough space to get a matrix blockrow. A row length of \n");
         printf("%d was not sufficient\n",(*values_space-1)/2);
         fflush(stdout);
         exit(1);
      }
      for (i = 0; i < index2; i++) t3[i] = (*values)[i];
      ML_free(*values);
      *values = t3;
   }

   /*iterate over the blocking data and store it*/
   bindx_old = &(input_matrix->bindx[bpntr_rows]);
   indx_old = &(input_matrix->indx[bpntr_rows]);
   if(index == 0){
     **indx = 0;
     /*indx_old++;*/
   }
   bindx_ptr = &((*bindx)[index]);
   indx_ptr = &((*indx)[index+1]);
   for(i = 0; i < *blocks; i++)
   {
     *bindx_ptr++ = *bindx_old++;
     *indx_ptr = indx_ptr[-1] - *indx_old++;
     *indx_ptr++ += *indx_old;
   }

   /*if we need to convert to global ids this will only be with the right B matrix*/
   if(place_holder->getrow->use_loc_glob_map == ML_YES)
   {
     col_size = input_matrix->cpntr[1]-input_matrix->cpntr[0];
     for(i = -*blocks; i < 0; i++)
     {
       bindx_ptr[i] = place_holder->getrow->loc_glob_map[bindx_ptr[i]*col_size]/col_size;
     }
   }

   /*iterate over the data and copy it*/
   val = &(input_matrix->val[start]);
   values_ptr = &((*values)[index2]);
   for(i = start; i < finish; i++)
   {
     *values_ptr++ = *val++;
   }
   return(1);
}
/*********************************************************************/
/* Get some matrix rows ( requested_rows[0 ... N_requested_rows-1] ) */
/* from the user's matrix and return this information  in            */
/* 'row_lengths, columns, values'.                                   */
/* If there is not enough space to complete this operation, return 0.*/
/* Otherwise, return 1.                                              */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/* data             On input, points to user's data containing       */
/*                  matrix values.                                   */
/* N_requested_rows On input, number of rows for which nonzero are   */
/*                  to be returned.                                  */
/* requested_rows   On input, requested_rows[0...N_requested_rows-1] */
/*                  give the row indices of the rows for which       */
/*                  nonzero values are returned.                     */
/* row_lengths      On output, row_lengths[i] is the number of       */
/*                  nonzeros in the row 'requested_rows[i]'          */
/*                  ( 0 <= i < N_requested_rows). NOTE: this         */
/*                  array is of size 'N_requested_rows'.             */
/* columns,values   On output, columns[k] and values[k] contains the */
/*                  column number and value of a matrix nonzero where*/
/*                  all the nonzeros for requested_rows[0] appear    */
/*                  first followed by the nonzeros for               */
/*                  requested_rows[1], etc. NOTE: these arrays are   */
/*                  of size 'allocated_space'.                       */
/* allocated_space  On input, indicates the space available in       */
/*                  'columns' and 'values' for storing nonzeros. If  */
/*                  more space is needed, return 0.                  */
/*********************************************************************/

int MSR_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   int    *bindx, j, row, start, finish;
   double *val;
   struct ML_CSR_MSRdata *input_matrix;
   ML_Operator *mat_in;

   mat_in = (ML_Operator *) data;
  input_matrix = (struct ML_CSR_MSRdata *) ML_Get_MyGetrowData(mat_in);

  bindx  = input_matrix->columns;
  val    = input_matrix->values;

  row    = *requested_rows;
  start  = bindx[row];
  finish = bindx[row+1];
  *row_lengths = finish - start + 1;
  if (*row_lengths > allocated_space) {
    ML_avoid_unused_param( (void *) &N_requested_rows);
    return(0);
  }

  /* diagonal */

  *columns++ = row;
  *values++  = val[row];

  /* off-diagonals */

  bindx = &(bindx[start]);
  val   = &(val[start]);

  for (j = start; j < finish; j++) {
     *columns++ = *bindx++;
     *values++  = *val++;
  }
   return(1);
}

int MSR_get_ones_rows(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   int    *bindx, j, row, start, finish;
   struct ML_CSR_MSRdata *input_matrix;
   ML_Operator *mat_in;

   mat_in = (ML_Operator *) data;
  input_matrix = (struct ML_CSR_MSRdata *) ML_Get_MyGetrowData(mat_in);
  bindx  = input_matrix->columns;

  row    = *requested_rows;
  start  = bindx[row];
  finish = bindx[row+1];
  *row_lengths = finish - start + 1;
  if (*row_lengths > allocated_space) {
    ML_avoid_unused_param( (void *) &N_requested_rows);
    return(0);
  }

  /* diagonal */

  *columns++ = row;
  *values++  = 1.;

  /* off-diagonals */

  bindx = &(bindx[start]);

  for (j = start; j < finish; j++) {
     *columns++ = *bindx++;
     *values++  = 1.;
  }
   return(1);
}

/*gets a point row from a VBR matrix same as above header this makes sure the matrix can go straight to a csr when converted.  Zeros are left in to allow for matrices already in VBR format to be put back into VBR more easily after they are passed to other processors as submatrices.*/
int VBR_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct ML_vbrdata *input_matrix;
   int  i, j, offset, point_rows, iminus1;
   int *cpntr, *bindx, *rpntr, *indx, *bpntr, startblock, endblock;
   double *val, cur_val;


   /*set alliases*/
   input_matrix = (struct ML_vbrdata *) ML_Get_MyGetrowData(data);
   bindx  = input_matrix->bindx;
   val    = input_matrix->val;
   rpntr  = input_matrix->rpntr;
   cpntr  = input_matrix->cpntr;
   bpntr  = input_matrix->bpntr;
   indx   = input_matrix->indx;

   /*find out where the requested row within the block structure.*/
   for(i = 0; rpntr[i] <= *requested_rows; i++);
   iminus1 = i-1;
   offset = indx[bpntr[iminus1]] + *requested_rows - rpntr[iminus1];
   point_rows = rpntr[i] - rpntr[iminus1];

   /*Make sure we have enough allocated space here*/
   if ((indx[bpntr[i]] - indx[bpntr[iminus1]])/point_rows > allocated_space) {
     ML_avoid_unused_param( (void *) &N_requested_rows);
     return(0);
   }

   /*Which blocks have data in them and how many of them there are*/
   startblock = bpntr[iminus1];
   endblock = bpntr[i];

   *row_lengths = 0; /*data points stored*/

   /*we need something to check here but in time we'll have it*/
   if(data->getrow->columns_loc_glob == ML_GLOBAL_INDICES)
   {
     /*iterate over the blocks*/
     for(i = startblock; i < endblock; i++)
     {
       /*iterate over columns in each block since the cuts are all the same just
         use the first pointer values*/
       for(j = cpntr[0]; j < cpntr[1]; j++)
       {
         cur_val = val[point_rows*row_lengths[0]+offset];
         values[*row_lengths] = cur_val;
         columns[*row_lengths] = j+bindx[i]*(cpntr[1]-cpntr[0]);
         row_lengths[0]++;
       }
     }
   }
   else
   {
     /*iterate over the blocks*/
     for(i = startblock; i < endblock; i++)
     {
       /*iterate over columns in each block*/
       for(j = cpntr[bindx[i]]; j < cpntr[bindx[i]+1]; j++)
       {
         cur_val = val[point_rows*row_lengths[0]+offset];
         values[*row_lengths] = cur_val;
         columns[*row_lengths] = j;
         row_lengths[0]++;
       }
     }
   }
   return(1);
}


int CSR_getrows(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
   register double *val;
   struct ML_CSR_MSRdata *in_mat;

   in_mat = (struct ML_CSR_MSRdata *) data;
   j = in_mat->rowptr[*requested_rows];
   *row_lengths = in_mat->rowptr[*requested_rows+1] - j;
   if (*row_lengths > allocated_space) return(0);

   bindx  = in_mat->columns + j;
   val    = in_mat->values + j;
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
      *values++  = *val++;
   }
   return(1);
}

/*gets one row of a CSR matrix ~IK*/
int CSR_getrow(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
   register double *val;
   struct ML_CSR_MSRdata *in_mat;

   in_mat = (struct ML_CSR_MSRdata *) data->data;
   j = in_mat->rowptr[*requested_rows];
   *row_lengths = in_mat->rowptr[*requested_rows+1] - j;
   if (*row_lengths > allocated_space) return(0);

   bindx  = in_mat->columns + j;
   val    = in_mat->values + j;
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
      *values++  = *val++;
   }
   return(1);
}

int CSR_get_one_row(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
   struct ML_CSR_MSRdata *in_mat;

   in_mat = (struct ML_CSR_MSRdata *) data->data;
   j = in_mat->rowptr[*requested_rows];
   *row_lengths = in_mat->rowptr[*requested_rows+1] - j;
   if (*row_lengths > allocated_space) return(0);

   bindx  = in_mat->columns + j;
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
      *values++  = 1.0;
   }
   return(1);
}

/* single precision version to save storage */
int sCSR_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
   register float *val;
   struct ML_CSR_MSRdata *inmat;

   inmat = (struct ML_CSR_MSRdata *) data->data;
   j = inmat->rowptr[*requested_rows];
   *row_lengths = inmat->rowptr[*requested_rows+1] - j;
   if (*row_lengths > allocated_space) return(0);

   bindx  = inmat->columns + j;
   val    = ((float *) inmat->values) + j;
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
      *values++  = (double) *val++;
   }
   return(1);
}

/* really cheap matrix format for null space within   */
/* edge element AMG (which contains only + or -1's).  */
int cCSR_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
   int     *rowptr,  row, itemp;
   register char    *val;
   struct ML_CSR_MSRdata *input_matrix;
   double sgn[3] = {0.,1.,-1.};
   ML_Operator *mat_in;

   mat_in = (ML_Operator *) data;
   row            = *requested_rows;
   input_matrix = (struct ML_CSR_MSRdata *) ML_Get_MyGetrowData(mat_in);
   rowptr = input_matrix->rowptr;
   if (rowptr != NULL) {
     itemp = rowptr[row];
     *row_lengths = rowptr[row+1] - itemp;
   }
   else {
     itemp = row;
     *row_lengths = 1;
   }
   if (*row_lengths > allocated_space) {
    ML_avoid_unused_param( (void *) &N_requested_rows);
    return(0);
  }

   bindx  = &(input_matrix->columns[itemp]);
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
   }
   val    = (char *) input_matrix->values;
   val    = &(val[itemp]);

   for (j = 0 ; j < *row_lengths; j++) {
     *values++ = sgn[ (int) *val++ ];
   }
   return(1);
}


/*********************************************************************/
/* Get some matrix rows ( requested_rows[0 ... N_requested_rows-1] ) */
/* from the user's matrix and return this information  in            */
/* 'row_lengths, columns, values'.  If there is not enough space to  */
/* complete this operation, return 0.  Otherwise, return 1.          */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/* data             On input, points to user's data containing       */
/*                  matrix values.                                   */
/* N_requested_rows On input, number of rows for which nonzero are   */
/*                  to be returned.                                  */
/* requested_rows   On input, requested_rows[0...N_requested_rows-1] */
/*                  give the row indices of the rows for which       */
/*                  nonzero values are returned.                     */
/* row_lengths      On output, row_lengths[i] is the number of       */
/*                  nonzeros in the row 'requested_rows[i]'          */
/*                  ( 0 <= i < N_requested_rows). NOTE: this         */
/*                  array is of size 'N_requested_rows'.             */
/* columns,values   On output, columns[k] and values[k] contains the */
/*                  column number and value of a matrix nonzero where*/
/*                  all the nonzeros for requested_rows[0] appear    */
/*                  first followed by the nonzeros for               */
/*                  requested_rows[1], etc. NOTE: these arrays are   */
/*                  of size 'allocated_space'.                       */
/* allocated_space  On input, indicates the space available in       */
/*                  'columns' and 'values' for storing nonzeros. If  */
/*                  more space is needed, return 0.                  */
/* getrows() function for a vector (i.e. matrix with 1 column).      */
/*********************************************************************/

int VECTOR_getrows(ML_Operator *mat_in, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   double *temp;
   int    i;
   ML_Operator *mat;

   mat = (ML_Operator *) mat_in;
   temp = (double *) ML_Get_MyGetrowData(mat);

   if (allocated_space < N_requested_rows) return(0);

   for (i = 0; i < N_requested_rows; i++) {
      row_lengths[i] = 1;
      columns[i]     = 0;
      values[i]      = temp[requested_rows[i]];
   }
   return(1);
}

/*********************************************************************/
/* Get some matrix rows ( requested_rows[0 ... N_requested_rows-1] ) */
/* from the user's matrix and return this information  in            */
/* 'row_lengths, columns, values'.  If there is not enough space to  */
/* complete this operation, return 0.  Otherwise, return 1.          */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/* data             On input, points to user's data containing       */
/*                  matrix values.                                   */
/* N_requested_rows On input, number of rows for which nonzero are   */
/*                  to be returned.                                  */
/* requested_rows   On input, requested_rows[0...N_requested_rows-1] */
/*                  give the row indices of the rows for which       */
/*                  nonzero values are returned.                     */
/* row_lengths      On output, row_lengths[i] is the number of       */
/*                  nonzeros in the row 'requested_rows[i]'          */
/*                  ( 0 <= i < N_requested_rows). NOTE: this         */
/*                  array is of size 'N_requested_rows'.             */
/* columns,values   On output, columns[k] and values[k] contains the */
/*                  column number and value of a matrix nonzero where*/
/*                  all the nonzeros for requested_rows[0] appear    */
/*                  first followed by the nonzeros for               */
/*                  requested_rows[1], etc. NOTE: these arrays are   */
/*                  of size 'allocated_space'.                       */
/* allocated_space  On input, indicates the space available in       */
/*********************************************************************/

int VBR_cnst_blk_getrows(ML_Operator *Amat_in, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct ML_vbrdata *input_matrix;
   int blk_row, count, N_rows, N_cols, row_offset, i, j, ii, start, k;
   int *cpntr, *bindx, *rpntr, *indx, *bpntr, step, row;
   double *val;


   input_matrix = (struct ML_vbrdata *) ML_Get_MyGetrowData(Amat_in);
   bindx  = input_matrix->bindx;
   val    = input_matrix->val;
   rpntr  = input_matrix->rpntr;
   cpntr  = input_matrix->cpntr;
   bpntr  = input_matrix->bpntr;
   indx   = input_matrix->indx;
   step = rpntr[1] - rpntr[0];


   count = 0;
   for (k = 0; k < N_requested_rows; k++) {
      row = requested_rows[k];
      blk_row = row/step;
      N_rows = rpntr[blk_row+1] - rpntr[blk_row];
      row_offset = row - rpntr[blk_row];
      row_lengths[k] = step*(bpntr[blk_row+1]-bpntr[blk_row]);
      if (count+row_lengths[k] > allocated_space) return(0);

      for (i = bpntr[blk_row]; i < bpntr[blk_row+1]; i++) {
         ii = bindx[i];
         N_cols = cpntr[ii+1] - cpntr[ii];
         start  = indx[i] + row_offset;

         for (j = 0; j < N_cols; j++) {
            columns[count] = j + cpntr[ii];
            values[count++] = val[start];
            start += N_rows;
         }
      }
   }
   return(1);
}

/*********************************************************************/
/* matvec in MSR format                                              */
/*********************************************************************/

int MSR_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{
  int i, j, Nrows, *bindx;
   double            *p2, *val, sum;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   ML_Comm           *comm;
   int *bindx_ptr;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat  = (ML_Operator *) Amat_in;
   comm  = Amat->comm;
   Nrows = Amat->matvec->Nrows;
   if ( ilen != olen && ilen != Nrows ) {
      printf("MSR_matvec error : lengths not matched.\n");
      exit(1);
   }
   temp  = (struct ML_CSR_MSRdata *) Amat->data;
   val   = temp->values;
   bindx = temp->columns;
   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      p2 = (double *) ML_allocate((Nrows+getrow_comm->minimum_vec_size+1)*
                                  sizeof(double));
      if (p2 == NULL)
	pr_error("MSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);

      for (i = 0; i < Nrows; i++) p2[i] = p[i];
      ML_exchange_bdry(p2,getrow_comm, Nrows, comm, ML_OVERWRITE,NULL);
   }
   else p2 = p;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
  j = bindx[0];
  bindx_ptr = &bindx[j];
  for (i = 0; i < Nrows; i++) {
    sum =  val[i]*p2[i];
    while (j+10 < bindx[i+1]) {
      sum += val[j+9]*p2[bindx_ptr[9]] +
	val[j+8]*p2[bindx_ptr[8]] +
	val[j+7]*p2[bindx_ptr[7]] +
	val[j+6]*p2[bindx_ptr[6]] +
	val[j+5]*p2[bindx_ptr[5]] +
	val[j+4]*p2[bindx_ptr[4]] +
	val[j+3]*p2[bindx_ptr[3]] +
	val[j+2]*p2[bindx_ptr[2]] +
	val[j+1]*p2[bindx_ptr[1]] +
	val[j]*p2[*bindx_ptr];
      bindx_ptr += 10;
      j += 10;
    }
    while (j < bindx[i+1]) {
      sum += val[j++] * p2[*bindx_ptr++];
    }
    ap[i] = sum;
  }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif


  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++) p[i] = p2[i];
     ML_free(p2);
  }
  return(1);
}

/*********************************************************************
* Matvec Av=y in CSR format.                                         *
*                                                                    *
* Amat_in   ML Operator (cast to void)                               *
* ilen      on input, length of x vector (int)                       *
* p         on input, x vector (int)                                 *
* olen      on input, length of y vector (int)                       *
* ap        on ouput, product Ax (double)                            *
*                                                                    *
*********************************************************************/

int CSR_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2, *val, sum, *ap2;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   /* Nrows   = Amat->outvec_leng; */
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("CSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);
     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("CSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);
   }
   else ap2 = ap;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   /* time only the real apply */
   t0 = GetClock();
#endif
   for (i = 0; i < Nstored; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
        sum  += val[k] * p2[bindx[k]];
     }

     ap2[i] = sum;
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

int sCSR_trans_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2/*, sum*/, *ap2;
   float *val;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = (float *) temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

      for (k = 0; k < i; k++) ap2[k] = 0.;
   }
   else {
     ap2 = ap;
     for (i = 0; i < olen; i++) ap2[i] = 0.;
   }


#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < ilen; i++) {
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
        ap2[bindx[k]] += ((double) val[k]) * p2[i];
     }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

int cCSR_trans_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2/*, sum*/, *ap2;
   char *val;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
   double            sgn[3] = {0.,1.,-1};
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = (char *) temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("cCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("cCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

      for (k = 0; k < i; k++) ap2[k] = 0.;
   }
   else {
     ap2 = ap;
     for (i = 0; i < olen; i++) ap2[i] = 0.;
   }

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   if (row_ptr != NULL) {
     for (i = 0; i < ilen; i++) {
       for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
	 {
	   ap2[bindx[k]] += (sgn[(int) val[k]]) * p2[i];
	 }
     }
   }
   else {
     for (i = 0; i < ilen; i++) {
       for (k = i; k < i+1; k++)
	 {
	   ap2[bindx[k]] += (sgn[(int) val[k]]) * p2[i];
	 }
     }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

int CSR_trans_ones_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2/*, sum*/, *ap2;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

      for (k = 0; k < i; k++) ap2[k] = 0.;
   }
   else {
     ap2 = ap;
     for (i = 0; i < olen; i++) ap2[i] = 0.;
   }


#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < ilen; i++) {
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
        ap2[bindx[k]] += 1.0 * p2[i];
     }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

int CSR_trans_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2/*, sum*/, *ap2;
   double *val;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = (double *) temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("sCSR_trans_matvec(%d): out of space\n",Amat->comm->ML_mypid);

      for (k = 0; k < i; k++) ap2[k] = 0.;
   }
   else {
     ap2 = ap;
     for (i = 0; i < olen; i++) ap2[i] = 0.;
   }


#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < ilen; i++) {
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
        ap2[bindx[k]] += ((double) val[k]) * p2[i];
     }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

int sCSR_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2, sum, *ap2;
   float             *val;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   /* Nrows   = Amat->outvec_leng; */
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = (float *) temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("sCSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("sCSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);

   }
   else ap2 = ap;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < Nstored; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
        sum  += ((double) val[k]) * p2[bindx[k]];
     }

     ap2[i] = sum;
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

/* really cheap matvec for null space within edge element AMG */
int cCSR_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2, sum, *ap2;
   char             *val;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
   double            sgn[3] = {0.,1.,-1.};
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   /* Nrows   = Amat->outvec_leng; */
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   val     = (char *) temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("cCSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("cCSR_matvec(%d): out of space\n",Amat->comm->ML_mypid);

   }
   else ap2 = ap;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   if (row_ptr != NULL) {
     for (i = 0; i < Nstored; i++) {
       sum = 0;
       for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
	 {
	   sum  += ( sgn[(int) val[ k]]) * p2[bindx[k]];
	 }

       ap2[i] = sum;
     }
   }
   else {
     for (i = 0; i < Nstored; i++) {
       sum = 0;
       for (k = i; k < i+1; k++)
	 {
	   sum  += ( sgn[(int) val[ k]]) * p2[bindx[k]];
	 }

       ap2[i] = sum;
     }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

/******************************************************************************/

int CSR_densematvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

           int i, jj, k, k2; /* Nrows,*/
           double            *p2, *val, sum, *ap2, *oldp2;
           struct ML_CSR_MSRdata *temp;
           ML_CommInfoOP     *getrow_comm;
           ML_Operator       *Amat;
           int               *row_ptr, Nstored;
           ML_Comm           *comm;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
           double t0;
#endif

           Amat    = (ML_Operator *) Amat_in;
           comm    = Amat->comm;
           /* Nrows   = Amat->outvec_leng; */
           Nstored = Amat->getrow->Nrows;
           temp    = (struct ML_CSR_MSRdata *) Amat->data;
           val     = temp->values;
           row_ptr = temp->rowptr;

           getrow_comm= Amat->getrow->pre_comm;
           if (getrow_comm != NULL) {
             p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                          sizeof(double));
             if (p2 == NULL)
               pr_error("CSR_dense_matvec(%d): out of space\n",Amat->comm->ML_mypid);

             for (i = 0; i < ilen; i++) p2[i] = p[i];

             ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

           }
           else p2 = p;

           getrow_comm= Amat->getrow->post_comm;
           if (getrow_comm != NULL) {
              i = Nstored+getrow_comm->minimum_vec_size + 1;
              if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
              ap2 = (double *) ML_allocate(i* sizeof(double));
              if (ap2 == NULL)
            pr_error("CSR_dense_matvec(%d): out of space\n",Amat->comm->ML_mypid);

           }
           else ap2 = ap;

           /*   jj = Amat->invec_leng; */

           oldp2 = p2;

           /* length of the product is the same for all vectors */
           /* It is enough to get the first length */
           k2    = row_ptr[1];

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
           t0 = GetClock();
#endif
           for (i = 0; i < Nstored; i++) {
             p2 = oldp2;
             sum = 0.;
             /* jj = 0 */
             for (k = 0; k < k2; k++) {
               /*       sum += val[j++]*p2[k];  */
               /*  sum += val[j++]*p2[jj++];  */

               sum += *val++ * *p2++;
             }
             /*  ap2[i] = sum;  */

             *ap2++ = sum;
           }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
           Amat->apply_without_comm_time += (GetClock() - t0);
#endif

           if (Amat->getrow->pre_comm != NULL) ML_free(oldp2);

           if (getrow_comm != NULL) {
              if (getrow_comm->remap != NULL) {
                 if (getrow_comm->remap_max != olen-1) {
                    printf("Error: The largest remapping index after communication\n");
                    printf("       should be one less than the vector's output\n");
                    printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
                    exit(1);
                 }
              }
              ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
              for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
              ML_free(ap2);
          }
          return(1);
}

/******************************************************************************/

int CSR_ones_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2, sum, *ap2;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
   double            *val_ptr;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   /* Nrows   = Amat->outvec_leng; */
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;
   val_ptr = temp->values;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("CSR_ones_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("CSR_ones_matvec(%d): out of space\n",Amat->comm->ML_mypid);

   }
   else ap2 = ap;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < Nstored; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
       if (val_ptr[k] != 0.0)
         sum  +=  p2[bindx[k]];
     }

     ap2[i] = sum;
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

/******************************************************************************/

int sCSR_ones_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, jj, k, /* Nrows,*/ *bindx;
   double            *p2, sum, *ap2;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   int               *row_ptr, Nstored;
   ML_Comm           *comm;
   float             *val_ptr;
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
#endif

   Amat    = (ML_Operator *) Amat_in;
   comm    = Amat->comm;
   /* Nrows   = Amat->outvec_leng; */
   Nstored = Amat->getrow->Nrows;
   temp    = (struct ML_CSR_MSRdata *) Amat->data;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;
   val_ptr = (float *) temp->values;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
     p2 = (double *) ML_allocate((getrow_comm->minimum_vec_size+ilen+1)*
                                  sizeof(double));
     if (p2 == NULL)
       pr_error("CSR_ones_matvec(%d): out of space\n",Amat->comm->ML_mypid);

     for (i = 0; i < ilen; i++) p2[i] = p[i];

     ML_exchange_bdry(p2,getrow_comm, ilen, comm, ML_OVERWRITE,NULL);

   }
   else p2 = p;

   getrow_comm= Amat->getrow->post_comm;
   if (getrow_comm != NULL) {
      i = Nstored+getrow_comm->minimum_vec_size + 1;
      if (getrow_comm->remap_max+1 > i) i = getrow_comm->remap_max+1;
      ap2 = (double *) ML_allocate(i* sizeof(double));
      if (ap2 == NULL)
	pr_error("CSR_ones_matvec(%d): out of space\n",Amat->comm->ML_mypid);

   }
   else ap2 = ap;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   t0 = GetClock();
#endif
   for (i = 0; i < Nstored; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
     {
       if (val_ptr[k] != 0.0)
         sum  +=  p2[bindx[k]];
     }

     ap2[i] = sum;
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   Amat->apply_without_comm_time += (GetClock() - t0);
#endif

   if (Amat->getrow->pre_comm != NULL) ML_free(p2);

   if (getrow_comm != NULL) {
      if (getrow_comm->remap != NULL) {
         if (getrow_comm->remap_max != olen-1) {
            printf("Error: The largest remapping index after communication\n");
            printf("       should be one less than the vector's output\n");
            printf("       length (%d vs %d)???\n",getrow_comm->remap_max,olen);
            exit(1);
         }
      }
      ML_exchange_bdry(ap2,getrow_comm, Nstored, comm, ML_ADD,NULL);
      for (jj = 0; jj < olen; jj++) ap[jj] = ap2[jj];
      ML_free(ap2);
  }
  return(1);
}

/******************************************************************************/

int localCSR_matvec(void *Amat_in, int ilen, double p[], int olen, double ap[])
{

   int i, k, /* Nrows,*/ *bindx;
   double            *p2, *val, sum, *ap2;
   struct ML_CSR_MSRdata *temp;
   int               *row_ptr, Nstored;

   Nstored = olen;
   temp    = (struct ML_CSR_MSRdata *) Amat_in;
   val     = temp->values;
   bindx   = temp->columns;
   row_ptr = temp->rowptr;

   p2 = p;
   ap2 = ap;
   if (ilen == -57) {
     ML_avoid_unused_param( (void *) &ilen);
   }


   for (i = 0; i < Nstored; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
        sum  += val[k] * p2[bindx[k]];

     ap2[i] = sum;
   }
  return(1);
}
/* ******************************************************************** */
/* constructor for DCSR matrix                                          */
/* -------------------------------------------------------------------- */

int ML_Matrix_DCSR_Create( ML_Matrix_DCSR **mat)
{
   ML_memory_alloc((void**) mat, sizeof(ML_Matrix_DCSR), "MA1");
   (*mat)->ML_id = ML_ID_MATRIXDCSR;
   (*mat)->mat_n = 0;
   (*mat)->mat_ia = NULL;
   (*mat)->mat_ja = NULL;
   (*mat)->mat_a  = NULL;
   (*mat)->comminfo = NULL;
   return 0;
}

/* ******************************************************************** */
/* destructor for DCSR matrix                                           */
/* -------------------------------------------------------------------- */

void ML_Matrix_DCSR_Destroy( ML_Matrix_DCSR *mat )
{
   if ( (mat)->mat_ia != NULL ) ML_free((mat)->mat_ia);
   if ( (mat)->mat_ja != NULL ) ML_free((mat)->mat_ja);
   if ( (mat)->mat_a  != NULL ) ML_free((mat)->mat_a);
   if ( (mat)->comminfo != NULL )
      ML_CommInfoOP_Destroy( &((mat)->comminfo) );
}

/* ******************************************************************** */
/* set data field                                                       */
/* -------------------------------------------------------------------- */

int ML_Matrix_DCSR_Set( ML_Matrix_DCSR *mat, int leng, int *ia, int *ja,
                        double *val)
{
   mat->mat_n = leng;
   mat->mat_ia = ia;
   mat->mat_ja = ja;
   mat->mat_a  = val;
   return 0;
}

/* ******************************************************************** */
/* set communication information field                                  */
/* -------------------------------------------------------------------- */

int ML_Matrix_DCSR_Set_Comm(ML_Matrix_DCSR *mat,ML_CommInfoOP *comminfo,
                            ML_Comm *comm)
{
   mat->comminfo = comminfo;
   mat->comm     = comm;
   return 0;
}

/* ********************************************************************** */
/* get a single row                                                       */
/* ---------------------------------------------------------------------- */

int ML_Matrix_DCSR_Getrow(ML_Operator *mat_in, int N_req_rows, int req_rows[],
        int allocated, int columns[], double values[], int row_lengths[])
{
   int             *rowptr,  row, itemp;
   register int    i, *bindx;
   register double *val;
   ML_Matrix_DCSR  *Amat;
   ML_Operator     *mat;

   mat = (ML_Operator *) mat_in;
   row  = *req_rows;
   Amat = (ML_Matrix_DCSR *) ML_Get_MyGetrowData(mat);
   rowptr = Amat->mat_ia;
   itemp  = rowptr[row];
   *row_lengths = rowptr[row+1] - itemp;
   if (*row_lengths > allocated) {
     ML_avoid_unused_param( (void *) &N_req_rows);
     return(0);
   }

   bindx  = &(Amat->mat_ja[itemp]);
   for (i = 0 ; i < *row_lengths; i++) *columns++ = *bindx++;
   val    = &(Amat->mat_a[itemp]);
   for (i = 0 ; i < *row_lengths; i++) *values++  = *val++;
   return(1);
}


/* ********************************************************************** */
/* matvec in CSR format                                                   */
/* ---------------------------------------------------------------------- */

int ML_Matrix_DCSR_Matvec(ML_Operator *mat_in,int ilen,double *x,int olen,double y[])
{
   int             i, k, Nrows, *bindx, nbytes;
   double          *y2, *val, sum;
   ML_CommInfoOP   *getrow_comm;
   int             *row_ptr;
   ML_Comm         *comm;
   ML_Matrix_DCSR  *Amat;
   ML_Operator     *mat;

   mat         = (ML_Operator *) mat_in;
   Amat        = (ML_Matrix_DCSR *) ML_Get_MyMatvecData(mat);

   comm        = Amat->comm;
   Nrows       = Amat->mat_n;
   val         = Amat->mat_a;
   bindx       = Amat->mat_ja;
   row_ptr     = Amat->mat_ia;
   getrow_comm = Amat->comminfo;

   if (olen != -57) {
     ML_avoid_unused_param( (void *) &olen);
   }

   if (getrow_comm != NULL)
   {
      nbytes = (getrow_comm->minimum_vec_size+ilen+1) * sizeof(double);
      y2 = (double *) ML_allocate( nbytes );
      for (i = 0; i < ilen; i++) y2[i] = x[i];
      ML_exchange_bdry(y2, getrow_comm, ilen, comm, ML_OVERWRITE,NULL);
   }
   else y2 = x;

   for (i = 0; i < Nrows; i++)
   {
      sum = 0;
      for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
         sum  += val[k] * y2[bindx[k]];
      y[i] = sum;
   }

   if (getrow_comm != NULL) ML_free(y2);

   return(1);
}

/************************************************************************/
/* Convert the data in csr_data from an MSR matrix to a CSR matrix.     */
/* Also, return the largest column number encountered.                  */
/*----------------------------------------------------------------------*/

int ML_MSR2CSR(struct ML_CSR_MSRdata *csr_data, int Nrows, int *Ncolumns)
{
  int  *row_ptr, *Tmat_bindx, i, j, lower, upper, nz_ptr, Ncols;
  double *Tmat_val, *diag;

  row_ptr = (int *) ML_allocate((Nrows+1)*sizeof(int));
  csr_data->rowptr  = row_ptr;
  Tmat_bindx = csr_data->columns;
  Tmat_val   = csr_data->values;

  diag    = (double *) ML_allocate(Nrows*sizeof(double));
  for (i = 0; i <= Nrows; i++) row_ptr[i] = Tmat_bindx[i];
  for (i = 0; i < Nrows; i++) diag[i] = Tmat_val[i];
  lower = row_ptr[0];
  row_ptr[0] = 0;
  nz_ptr = 0;
  Ncols = -1;
  for (i = 0; i < Nrows; i++) {

    upper = row_ptr[i+1];
    if ( diag[i] != 0.0) {
      Tmat_bindx[nz_ptr] = i;
      Tmat_val[nz_ptr++] = diag[i];
      if (i > Ncols) Ncols = i;
    }
    for (j = lower; j < upper; j++) {
      if (Tmat_val[j] != 0.0) {
    Tmat_bindx[nz_ptr] = Tmat_bindx[j];
    Tmat_val[nz_ptr++] = Tmat_val[j];
        if (Tmat_bindx[j] > Ncols) Ncols = Tmat_bindx[j];
      }
    }
    row_ptr[i+1] = nz_ptr;
    lower = upper;
  }
  Ncols++;
  ML_free(diag);
  *Ncolumns = Ncols;
  return 0;
}


#ifdef WKC
/* WKC -- note that the double * coming in are really Epetra_MultiVectors */

/*********************************************************************/
/* matvec in MSR format                                              */
/*********************************************************************/
int MSR_matvec_WKC(ML_Operator *Amat_in, int ilen, double *ep_p, int olen, double *ep_ap)
{
   int i, j, Nrows, *bindx;
   double            **big_p;
   double            *p2, *val, sum;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   ML_Comm           *comm;
   int *bindx_ptr;

   Amat  = (ML_Operator *) Amat_in;
   comm  = Amat->comm;
   Nrows = Amat->matvec->Nrows;
   if ( ilen != olen && ilen != Nrows ) {
      pr_error("MSR_matvec error : lengths not matched.\n");
   }
   temp  = (struct ML_CSR_MSRdata *) Amat->data;
   val   = temp->values;
   bindx = temp->columns;

   Epetra_MultiVector &X ( *(Epetra_MultiVector *)ep_p );
   Epetra_MultiVector &Y ( *(Epetra_MultiVector *)ep_ap );
   double **pp_p;
   double **pp_ap;
   X.ExtractView ( &pp_p );
   Y.ExtractView ( &pp_ap );

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      big_p = new double *[X.NumVectors()];
   for ( int KK = 0 ; KK != X.NumVectors() ; KK++ ) {
      big_p[KK] = (double *) ML_allocate((Nrows+getrow_comm->minimum_vec_size+1)*
                                  sizeof(double));
      for (i = 0; i < Nrows; i++) big_p[KK][i] = pp_p[KK][i];
      ML_exchange_bdry(big_p[KK],getrow_comm, Nrows, comm, ML_OVERWRITE,NULL);
   }
   }
   else {
    big_p = pp_p;
   }

  for ( int KK = 0 ; KK != X.NumVectors() ; KK++ ) {
  p2 = big_p[KK];
  double *ap = pp_ap[KK];

  j = bindx[0];
  bindx_ptr = &bindx[j];
  for (i = 0; i < Nrows; i++) {
    sum =  val[i]*p2[i];
    while (j+10 < bindx[i+1]) {
      sum += val[j+9]*p2[bindx_ptr[9]] +
        val[j+8]*p2[bindx_ptr[8]] +
        val[j+7]*p2[bindx_ptr[7]] +
        val[j+6]*p2[bindx_ptr[6]] +
        val[j+5]*p2[bindx_ptr[5]] +
        val[j+4]*p2[bindx_ptr[4]] +
        val[j+3]*p2[bindx_ptr[3]] +
        val[j+2]*p2[bindx_ptr[2]] +
        val[j+1]*p2[bindx_ptr[1]] +
        val[j]*p2[*bindx_ptr];
      bindx_ptr += 10;
      j += 10;
    }
    while (j < bindx[i+1]) {
      sum += val[j++] * p2[*bindx_ptr++];
    }
    ap[i] = sum;
  }
  }


  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++)
     for ( int KK = 0 ; KK != X.NumVectors() ; KK++ ) {
       X[KK][i] = big_p[KK][i];
     }
     for ( int KK = 0 ; KK != X.NumVectors() ; KK++ ) {
        ML_free ( big_p[KK] );
     }
     delete [] big_p;
  }
  return(1);
}

#endif
