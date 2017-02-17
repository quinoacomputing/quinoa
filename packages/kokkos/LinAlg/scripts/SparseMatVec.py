#@HEADER
# ************************************************************************
# 
#          Kokkos: Node API and Parallel Node Kernels
#              Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

'''Generate C++ code for sequential sparse matrix-(multi)vector multiply.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: June 2012

Introduction
============

This module generates C++ code for sequential sparse matrix-vector
multiply, where there may be one or more right-hand side vector(s).
We commonly abbreviate this as SpMV or SpMM (where the latter "M"
stands for "multiple vectors").  The module makes many routines, one
for each combination of parameters relating to the following:

- The sparse matrix format: compressed sparse row (CSR) or compressed
  sparse column (CSC).
- The dense matrix ("multivector") data layout: column or row major
- Whether the routine accepts a single input/output vector(s), or
  separate input and output vector(s).
- Whether to use the conjugate of each sparse matrix element.

Code variants
=============

Sparse and dense matrix storage formats
---------------------------------------

Both the CSC and CSR formats use the standard three-array 'ptr',
'ind', 'val' representation, where the 'ptr' array has one more entry
than the number of columns resp. rows.  All the routines take a dense
matrix ("multivector") with one or more columns as input, and
overwrite another dense matrix ("multivector") with the same number of
columns as output.  These dense matrices may be stored in either
column-major (Fortran style) or row-major (C style) order.  We assume
that ptr, ind, and val encode exactly the matrix to use; there are no
special cases for assuming an implicit unit diagonal or only using the
lower or upper triangle.  We do support these in the sparse triangular
solve case, but not in the sparse matrix-vector multiply case.

Conjugate of each sparse matrix element
---------------------------------------

This feature lets us implement CSR conjugate transpose by treating the
CSR three arrays as a CSC-format sparse matrix.

Properties of the generated C++ code
====================================

Expected arguments
------------------

The sparse triangular solve routines all call the output vector(s) X.
The out-of-place routines call the input vector(s) Y.  This follows
the analogy of sparse matrix-vector multiply: Y = A * X, so X = A \ Y
(using Matlab's backslash notation for "solve, like X = inv(A)*Y but
without computing the inverse of A").

Template parameters
-------------------

Each generated out-of-place routine has four template parameters:
Ordinal, MatrixScalar, DomainScalar, and RangeScalar.  Ordinal is the
type of indices to use, MatrixScalar the type of entries in the sparse
matrix, DomainScalar the type of entries in the input (dense) matrix,
and RangeScalar the type of entries in the output (dense) matrix.  The
in-place routines omit DomainScalar, since the input vector(s) is/are
also the output vector(s).

Output vector(s) do(es) not overwrite the input vector(s)
---------------------------------------------------------

CSR sparse matrix-vector multiply is best suited to "out-of-place"
computation, where the output vector and input vector are separate.
Since we favor CSR and include CSC only for the transpose case, we
choose only to generate out-of-place versions of CSC.

Expected performance
--------------------

We assume the following performance characteristics of sparse
matrix-vector multiply:

1. Reading the entries (indices and values) of the sparse matrix is
   the main cost.
2. Avoid branches whenever possible.

Hard-coding each routine to its set of options avoids branches in
inner loops, which should result in faster code.  The generated C++
code favors cache-based CPU architectures.  It assumes that the main
cost of sequential sparse matrix-vector multiply is reading the sparse
matrix.  Thus, all routines amortize the cost of reading the sparse
matrix over all the columns of the input and output matrices.  This
introduces little or no additional cost if there is only one column in
the input and output matrices.  For multiple columns, this should
always pay off for row-major storage.  For column-major storage, this
should pay off as long as the sparse matrix has on average more
entries per row than the number of MatrixScalar (see below) values
that can fit in a cache line.  Row-major storage should generally be
faster than column-major storage if there are multiple input and
output columns.  (This is because column-major storage accesses the
input and output matrices with nonunit stride.)

Note that sparse matrix data structures other than compressed sparse
row or column often perform much better, especially if you make
assumptions about the structure of the sparse matrix.  Rich Vuduc's
PhD thesis, the OSKI project, etc. all refer to this.  We do not
attempt to generate such code here.

Algorithm variant
-----------------

The usual CSC and CSR sparse matrix-vector multiply routines have two
nested 'for' loops: one for the columns resp. rows, and one for the
entries within a column resp. row.  We have chosen instead to generate
a single 'for' loop over all the entries in the sparse matrix.  Within
that 'for' loop is a 'while' loop for incrementing the current column
resp. row index.  Epetra uses this variant, and our experience is that
it performs slightly better.

Here is a sketch of a correctness proof for this variant (I'll
consider the CSR case without loss of generality):

Invariants inside the while loop:
* 0 <= k < ptr[numRows]
  
Invariants inside this loop, before ++i
* 0 <= i < numRows
* k >= ptr[i+1], which means that i is still too small.
* We have not yet initialized Y(i+1,:).
  
Since we know that 0 <= k < ptr[numRows], we know that A_ij = val[k]
and j = ind[k] are valid.  Thus, the correct i is the one for which
ptr[i] <= k < ptr[i+1].  If ptr[i] == ptr[i+1], then the corresponding
row i is empty.  In that case, k >= ptr[i] and k >= ptr[i+1] as well,
so this loop will move i past that row.
  
If the last row of the matrix is empty, then ptr[numRows-1] ==
ptr[numRows].  However, k < ptr[numRows] always (see above invariant),
so we would never enter this 'while' loop in that case.  Thus, we
don't need to check in the 'while' clause whether i < numRows.
  
We need a while loop, and not just an if test, specifically for the
case of empty rows.  If we forbid empty rows (this is easy to do by
simply adding an entry with a zero value to each empty row when
constructing ptr,ind,val), then we can replace the while loop with a
single if test.  This saves a branch.

How to use the code generator
=============================

Normal (nonexpert) use
----------------------

If you run this module as an executable script, like this:

$ python SparseMatVec.py

it will write two header files to the current directory.
Kokkos_Raw_SparseMatVec_decl.hpp will contain function declarations,
and Kokkos_Raw_SparseMatVec_def.hpp will contain function definitions
(see below).  Users who do not want to modify the generated code at
all or change the output file names or output directory will use this
script in that way.

Expert use
----------

All functions that generate code return the C++ code as a string.  The
returned string itself is a "raw" string in the Python sense.  If you
want to read it at the Python prompt, 'print' the string.

The top-level functions in this module, emitHeaderDeclFile() and
emitHeaderDefFile(), create entire header files with all possible code
variants.  The 'Decl' version makes a header file with just the
declarations of the functions, whereas the 'Def' version makes a
header file with their definitions.  We separate these so that it
won't be hard to implement explicit instantiation later, in case this
is useful.

If you want to make a single function declaration, including its
documentation, call emitFuncDecl().  For a single function definition
(without documentation, which belongs to the declaration in any case),
call emitFuncDef().  emitFuncDoc() takes the same arguments, and
generates Doxygen-formatted documentation for the routine that would
be generated by emitFuncDef() with the same input dictionary (see
below).

My intent is for the generated code to be dumped to a file, compiled
by a C++ compiler, and then linked into C++ program at link time,
before that program runs.  However, if you are adventurous, you might
like to try generating code at run time and patching it into a running
program.  (See the "Related work" section for details.)  I make no
promises about the suitability of this code for that purpose.

Parameters used to generate routines
------------------------------------

Many functions in this module take a dictionary 'defDict' as input.
The dictionary defines which variant of sparse triangular solve to
generate.  It must have the following fields:

sparseFormat: Sparse matrix storage format.  Currently supported
  formats: 'CSR' for compressed sparse row, 'CSC' for compressed
  sparse column.

dataLayout: Layout of the dense multivectors which are the input and
  output arguments of the sparse matrix-vector multiply routines.  The
  current values we accept are 'column major' or 'row major'.

conjugateMatrixEntries: Whether to compute with the (complex)
  conjugate of the matrix entries before using them.  We use this
  option to implement sparse matrix-vector multiply with the conjugate
  transpose of the matrix.  If the MatrixScalar type (the type of
  entries in the sparse matrix) is real, then the option does nothing.

It may also have the following fields:

hardCodeNumVecs (Boolean): Whether to hard-code the number of columns
  in the input and output multivectors to numVecs (see below).
  Default is False.

numVecs (nonnegative integer): If 'hardCodeNumVecs' (see above) is
  True, hard-code the number of columns in the input and output
  multivectors to this value.  Ignored if hardCodeNumVecs=False.

unrollLength (positive integer): For hardCodeNumVecs==False, this
  specifies the strip-mine length for unrolling updates over the input
  and output multivectors.  (This is ignored if hardCodeNumVecs==True,
  since in that case, we simply unroll over all columns.)  We've
  written the strip-mined loops so that they are correct for any
  nonnegative number of columns (numVecs) in the input and output
  multivectors.

Related work
============

Code that writes code is not a particularly new idea, even in sparse
matrix computations.  One of the earliest examples I know would take a
sparse matrix structure and generate a custom sparse factorization
code for all matrices with that structure.

The SEJITS (Selective Embedded Just-In-Time Specialization) effort at
the University of California Berkeley's Computer Science department
(point of contact: Prof. Armando Fox) lets programmers write in a
high-level language like Python or Ruby.  A "specializer" then
translates their code at run time to use generated code and/or
optimized routines in a lower-level language.  Our work differs from
theirs because we're only using Python to generate code; we don't
include source-to-source translation from Python into C++, and we
don't intend end consumers of the generated routines to call them from
Python.'''

from string import Template
from os.path import basename
from kokkos import makeCopyrightNotice
from SparseCodeGen import emitDenseAref, emitDenseArefFixedCol


def makeDefDict (sparseFormat, dataLayout, conjugateMatrixEntries, hardCodeNumVecs=False, numVecs=1, unrollLength=4):
    '''Make a suitable input dictionary for any function here that takes one.

    This function is mainly useful for interactive debugging.  See the
    "Parameters used to generate routines" section in this module's
    documentation for an explanation of this function's arguments.'''
    return {'sparseFormat': sparseFormat,
            'dataLayout': dataLayout,
            'conjugateMatrixEntries': conjugateMatrixEntries,
            'hardCodeNumVecs': hardCodeNumVecs,
            'numVecs': numVecs,
            'unrollLength': unrollLength}

def emitFuncDeclVariants (indent):
    '''Generate declarations of all sparse matrix-vector multiply variants.'''
    return emitFuncVariants (emitFuncDecl, indent)

def emitFuncDefVariants (indent):
    '''Generate definitions of all sparse matrix-vector multiply variants.'''
    return emitFuncVariants (emitFuncDef, indent)

def makesSense (defDict):
    '''Whether the sparse matrix-vector multiply variant specified by defDict makes sense.

    "Makes sense" means that the combination of parameters specified
    by the input dictionary results in a valid variant.  We don't
    generate variants for which this function returns False.'''
    return True

def emitFuncVariants (f, indent):
    '''Return a string with all sparse matrix-vector multiply variants.

    f: a function that takes (defDict, indent) and returns a string.
    
    indent: nonnegative integer indent level.  All lines of code
      emitted by this function are indented by this many spaces.'''

    maxHardCodedNumVecs = 4
    s = ''
    for conjugateMatrixEntries in [False, True]:
        for dataLayout in ['column major', 'row major']:
            for sparseFormat in ['CSC', 'CSR']:
                for hardCodeNumVecs in [False, True]:
                    if hardCodeNumVecs:
                        unrollLength = 4 # ignored in this case
                        for numVecs in xrange (1, maxHardCodedNumVecs+1):
                            d = makeDefDict (sparseFormat, dataLayout, conjugateMatrixEntries, hardCodeNumVecs, numVecs, unrollLength)
                            if makesSense (d):
                                s = s + f(d, indent) + '\n'
                    else:
                        numVecs = 1 # ignored in this case
                        for unrollLength in [1, 4]: # unrollLength==1 means don't unroll loops.
                            d = makeDefDict (sparseFormat, dataLayout, conjugateMatrixEntries, hardCodeNumVecs, numVecs, unrollLength)
                            if makesSense (d):
                                s = s + f(d, indent) + '\n'
    return s
    
def emitFuncName (defDict):
    '''Emit the function's name.'''
    # Model:
    # matVecCsrColMajorConj

    name = ''
    name = name + 'matVec'
    # Sparse matrix storage format.
    if defDict['sparseFormat'] == 'CSC':
        name = name + 'Csc'
    elif defDict['sparseFormat'] == 'CSR':    
        name = name + 'Csr'
    else:
        raise ValueError('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')
    # Layout of the dense input and output (multi)vectors.
    if defDict['dataLayout'] == 'column major':
        name = name + 'ColMajor'
    elif defDict['dataLayout'] == 'row major':    
        name = name + 'RowMajor'
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')
    # Various Boolean options
    if defDict['conjugateMatrixEntries']:
        name = name + 'Conj'
    if defDict['hardCodeNumVecs']:
        name = name + str (defDict['numVecs']) + 'Vec'
    elif defDict['unrollLength'] > 1: # unrollLength == 1 means we don't unroll loops
        name = name + str (defDict['unrollLength']) + 'Unrolled'
    return name

def emitFuncDecl (defDict, indent=0):
    '''Emit the function's declaration, including documentation.'''
    return emitFuncDoc(defDict, indent) + '\n' + emitFuncSig(defDict, indent) + ';\n'

def emitFuncDef (defDict, indent=0):
    '''Emit the function's definition.'''
    return emitFuncSig(defDict, indent) + '\n' + \
        emitFuncBody (defDict, indent)

def emitFuncSig (defDict, indent=0):
    '''Emit the function's signature (without terminal end-of-line).'''
    sig = ''
    ind = ' '*indent
    sig = sig + ind + 'template<class Ordinal,\n' + \
        ind + '         class MatrixScalar,\n' + \
        ind + '         class DomainScalar,\n' + \
        ind + '         class RangeScalar>\n' + \
        ind + 'void\n' + \
        ind + emitFuncName(defDict) + ' (\n' + \
        ind + '  const Ordinal numRows,\n' + \
        ind + '  const Ordinal numCols,\n'
    if not defDict['hardCodeNumVecs']:
        sig = sig + ind + '  const Ordinal numVecs,\n'
    sig = sig + \
        ind + '  const RangeScalar& beta,\n' + \
        ind + '  RangeScalar Y[],\n' + \
        ind + '  const Ordinal ${denseRowCol}StrideY,\n' + \
        ind + '  const RangeScalar alpha,\n' + \
        ind + '  const Ordinal ptr[],\n' + \
        ind + '  const Ordinal ind[],\n' + \
        ind + '  const MatrixScalar val[],\n' + \
        ind + '  const DomainScalar X[],\n' + \
        ind + '  const Ordinal ${denseRowCol}StrideX)'

    if defDict['dataLayout'] == 'column major':
        denseRowCol = 'col'
    elif defDict['dataLayout'] == 'row major':
        denseRowCol = 'row'        
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')
    return Template(sig).substitute(denseRowCol=denseRowCol)

def emitFuncBody (defDict, indent=0):
    '''Generate the sparse matrix-vector multiply function body, including { ... }.'''
    ind = ' '*indent
    body = ''
    body = body + \
        ind + '{\n' + \
        ind + '  ' + 'typedef Teuchos::ScalarTraits<RangeScalar> STS;\n\n'
    if defDict['hardCodeNumVecs']:
        body = body + ind + '  ' + 'const Ordinal numVecs = ' + str (defDict['numVecs']) +';\n'
    if defDict['sparseFormat'] == 'CSC':
        # CSC requires prescaling the output vector(s) Y.
        body = body + emitPreScaleLoop (defDict, indent+2)
        loopIndex = 'j'
        otherIndex = 'i'
    elif defDict['sparseFormat'] == 'CSR':
        loopIndex = 'i'
        otherIndex = 'j'
    else:
        raise ValueError ('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')

    return body + \
        emitForLoopPreface (defDict, loopIndex, indent+2) + \
        ind + ' '*2 + 'if (alpha == STS::zero()) {\n' + \
        ind + ' '*4 + 'return; // Our work is done!\n' + \
        ind + ' '*2 + '}\n' + \
        emitForLoop (defDict, indent+2) + \
        ' '*indent + '}\n'

def emitPreScaleLoop (defDict, indent=0):
    '''Emit the prescale loop for sparse matrix-vector multiply.

    Sparse matrix-vector multiply (SpMV) computes Y := beta*Y +
    alpha*A*X, for scalar constants alpha and beta.  Some
    implementations of SpMV need to start with Y := beta*Y.  We call
    this a "prescale."  This function generates the prescale loop.
    SpMV must pass over the entries of Y anyway, so prescaling is
    suboptimal with respect to the number of reads and writes of the
    entries of Y.  However, it's necessary sometimes.

    Prescaling has two special cases: beta = 0, and beta = 1.  If beta
    = 1, then a full prescale isn't necessary, because we already
    formulate the sparse matrix-vector product as an update (Y(i) =
    Y(i) + alpha * A(i,j) * X(j)).  If beta = 0, then we can simplify
    the prescale by replacing the entries of Y with zeros.  (We follow
    the Sparse BLAS convention that the result of a scale is zero if
    beta is zero, regardless of any NaN or Inf entries in the vector.)
    '''

    ind = ' ' * indent
    layout = defDict['dataLayout']
    if layout != 'column major' and layout != 'row major':
        raise ValueError ('Invalid dataLayout "' + layout + '"')

    s = ind + '// Prescale: Y := beta * Y.\n' + \
        ind + 'if (beta == STS::zero()) {\n'
    if layout == 'column major':
        Y_j = emitDenseAref (defDict, 'Y', '0', 'j')
        s = s + \
            ind + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n' + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*6 + 'Y_j[i] = STS::zero();\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        s = s + \
            ind + ' '*2 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
            ind + ' '*4 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*6 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*6 + 'Y_i[j] = STS::zero();\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    s = s + \
        ind + '}\n' + \
        ind + 'else if (beta != STS::one()) {\n'
    # It's more efficient to put stride-1 access in the inner loop.
    if layout == 'column major':
        Y_j = emitDenseAref (defDict, 'Y', '0', 'j')
        s = s + \
            ind + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n' + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + 'Y_j[i] = beta * Y_j[i];\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        s = s + \
            ind + ' '*2 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
            ind + ' '*4 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*6 + 'Y_i[j] = beta * Y_i[j];\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    return s + ind + '}\n'

def emitForLoopPreface (defDict, loopIndex, indent=0):
    ind = ' '*indent
    s = ind + 'Ordinal ' + loopIndex + ' = 0;\n'
    if defDict['sparseFormat'] == 'CSR':
        Y_0c = emitDenseAref (defDict, 'Y', '0', 'c')
        # No real need to unroll this loop, since it's only for the
        # first entry of Y.  We have to treat beta==0 separately in
        # order to follow the Sparse BLAS convention to replace Inf
        # and NaN entries in Y with zero if beta==0, rather than
        # allowing them to propagate according to IEEE 754.
        s = s + \
            ind + '// Special case for CSR only: Y(0,:) = 0.\n' + \
            ind + 'if (beta != STS::zero()) {\n' + \
            ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
            ind + ' '*4 + Y_0c + ' = beta * ' + Y_0c + ';\n' + \
            ind + ' '*2 + '}\n' + \
            ind + '}\n' + \
            ind + 'else {\n' + \
            ind + ' '*2 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
            ind + ' '*4 + Y_0c + ' = STS::zero();\n' + \
            ind + ' '*2 + '}\n' + \
            ind + '}\n'
    return s

def emitForLoopImpl (defDict, alphaIsOne, indent=0):
    '''Generate the sparse mat-vec 'for' loop, for a specific alpha case (alpha == 1 or alpha != 1).

    defDict: The usual dictionary.
    alphaIsOne (Boolean): True if alpha == 1, False if alpha != 1.
    indent (integer): How many spaces to indent each line of emitted code.'''

    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']

    X_jc = emitDenseAref (defDict, 'X', 'j', 'c')
    Y_ic = emitDenseAref (defDict, 'Y', 'i', 'c')

    if defDict['sparseFormat'] == 'CSC':
        loopIndex = 'j'
        otherIndex = 'i'
        RowCol = 'Col'
    else:
        loopIndex = 'i'
        otherIndex = 'j'
        RowCol = 'Row'

    if defDict['conjugateMatrixEntries']:
        getMatVal = 'Teuchos::ScalarTraits<MatrixScalar>::conjugate (val[k])'
    else:
        getMatVal = 'val[k]'

    ind = ' '*indent
    s = ''
    s = s + \
        ind + 'for (Ordinal k = 0; k < nnz; ++k) {\n' + \
        ind + ' '*2 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
        ind + ' '*2 + 'const Ordinal ${otherIndex} = ind[k];\n'
    # Here comes the 'while' loop for updating the row index (for CSR,
    # or column index for CSC).  This is where the one-for-loop
    # variant differs from the two-for-loop variant.
    s = s + \
        ind + ' '*2 + 'while (k >= ptr[${loopIndex}+1]) {\n' + \
        ind + ' '*4 + '++${loopIndex};\n'
    if defDict['sparseFormat'] == 'CSR':
        # CSR mat-vec can merge scaling Y by beta into the iteration
        # over rows.  CSC mat-vec can't do that; it has to prescale.
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        s = s + \
            ind + ' '*4 + '// We haven\'t seen row i before; scale Y(i,:) by beta.\n'
        if defDict['hardCodeNumVecs']:
            if defDict['numVecs'] > 1:
                s = s + ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
                for c in xrange (0, defDict['numVecs']):
                    Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                    s = s + ind + ' '*4 + Y_ic + ' *= beta;\n'
            else:
                Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
                s = s + ind + ' '*4 + Y_ic + ' *= beta;\n'
        elif unrollLength > 1:
            s = s + \
                ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
                ind + ' '*4 + 'Ordinal c = 0;\n' + \
                ind + ' '*4 + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
                ind + ' '*4 + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
                ind + ' '*4 + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'
                # Unrolled part of the loop.
            for c in xrange (0, unrollLength):
                Y_ic_fixed = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                s = s + ind + ' '*6 + Y_ic_fixed + ' *= beta;\n'
            Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
            # "Leftover" part of the loop.
            s = s + ind + ' '*4 + '}\n' + \
                ind + ' '*4 + 'for ( ; c < numVecs; ++c) {\n' + \
                ind + ' '*6 + Y_ic + ' *= beta;\n' + \
                ind + ' '*4 + '}\n'
        else: # unrollLength == 1, which means don't unroll loops at all.
            s = s + ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
            Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')            
            s = s + \
                ind + ' '*4 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                ind + ' '*6 + Y_ic + ' = beta * ' + Y_ic + ';\n' + \
                ind + ' '*4 + '}\n'
    s = s + ind + ' '*2 + '}\n' # End of the 'while' loop for advancing the row/column index.
    s = s + \
        emitUpdateLoop (defDict, alphaIsOne, indent+2) + \
        ind + '}\n'
    return Template(s).substitute (loopIndex=loopIndex, \
                                       otherIndex=otherIndex, \
                                       RowCol=RowCol, \
                                       getMatVal=getMatVal)

def emitForLoop (defDict, indent=0):
    '''Generate both branches of the sparse mat-vec 'for' loop, for each alpha case (alpha == 1 or alpha != 1).

    defDict: The usual dictionary.
    indent: The number of spaces to indent each line of emitted code; a nonnegative integer.'''

    if defDict['sparseFormat'] == 'CSC':
        RowCol = 'Col'
    else:
        RowCol = 'Row'

    ind = ' '*indent
    return ind + 'const Ordinal nnz = ptr[num' + RowCol + 's];\n' + \
        ind + 'if (alpha == STS::one()) {\n' + \
        emitForLoopImpl (defDict, True, indent+2) + \
        ind + '}\n' + \
        ind + 'else { // alpha != STS::one()\n' + \
        emitForLoopImpl (defDict, False, indent+2) + \
        ind + '}\n'

def emitUpdateLoop (defDict, alphaIsOne, indent=0):
    '''Return the update loop code for Y(i,:) = Y(i,:) + alpha*A_ij*X(j,:).

    alphaIsOne (Boolean): if True, then we know that alpha is one and
      don't have to multiply by alpha.  If False, we must multiply by
      alpha.
    indent: The number of spaces to indent each line of emitted code; a nonnegative integer.'''

    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']
    Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
    X_j = emitDenseAref (defDict, 'X', 'j', '0')

    ind = ' '*indent
    s = ''
    if defDict['hardCodeNumVecs']:
        s = s + emitUpdateLoopFixedNumVecs (defDict, alphaIsOne, indent)
    elif unrollLength > 1:
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        s = s + \
            ind + 'Ordinal c = 0;\n' + \
            ind + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
            ind + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
            ind + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'            
        # Unrolled part of the loop.
        for c in xrange (0, unrollLength):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if alphaIsOne:
                s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
            else:
                s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
        # "Leftover" part of the loop.
        X_jc = emitDenseAref (defDict, 'X_j', '0', 'c', strideName='X')
        Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
        s = s + ind + 'for ( ; c < numVecs; ++c) {\n'
        if alphaIsOne:
            s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
    else: # unrollLength == 1
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        X_jc = emitDenseAref (defDict, 'X_j', '0', 'c', strideName='X')
        Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
        s = s + ind + 'for (Ordinal c = 0; c < numVecs; ++c) {\n'
        if alphaIsOne:
            s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
    return s

def emitUpdateLoopFixedNumVecs (defDict, alphaIsOne, indent=0):
    '''Return fixed-numVecs update loop code for Y(i,:) += alpha * A_ij * X(j,:).

    The returned code completely unrolls the loop over all columns of Y.
    This is probably only a good idea for a small number of columns.

    defDict: The usual dictionary.
    alphaIsOne (Boolean): if True, then we know that alpha is one and
      don't have to multiply by alpha.  If False, we must multiply by
      alpha.
    indent: The number of spaces to indent; a nonnegative integer.'''
    if not defDict['hardCodeNumVecs']:
        raise ValueError('Only call this if defDict[\'hardCodeNumVecs\']==True.')
    numVecs = defDict['numVecs']
    X_j = emitDenseAref (defDict, 'X', 'j', '0', strideName='X')
    Y_i = emitDenseAref (defDict, 'Y', 'i', '0', strideName='Y')

    ind = ' '*indent
    s = ''

    if numVecs > 1:
        # Only make X_j and Y_i pointers if numVecs > 1.
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        for c in xrange(0, numVecs):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if alphaIsOne:
                line = ind + Y_ic + ' += A_ij * ' + X_jc + ';\n'
            else:
                line = ind + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
            s = s + line
    else: # numVecs == 1
        X_jc = emitDenseArefFixedCol (defDict, 'X', 'j', 0, strideName='X')
        Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
        if alphaIsOne:
            line = ind + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            line = ind + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + line
    return s


def emitFuncDoc (defDict, indent=0):
    '''Emit the sparse matrix-vector multiply routine's documentation.
    
    This generates the documentation (in Doxygen-compatible format)
    for a sparse matrix-vector multiply routine.'''

    sparseFormat = defDict['sparseFormat']
    dataLayout = defDict['dataLayout']
    conjugateMatrixEntries = defDict['conjugateMatrixEntries']

    if sparseFormat == 'CSC':
        fmtColRow = 'row'
        fmtRowCol = 'column'
        startIndex = 'startCol'
        endIndex = 'endColPlusOne'
        numIndices = 'numCols'
    elif sparseFormat == 'CSR':
        fmtColRow = 'column'
        fmtRowCol = 'row'
        startIndex = 'startRow'
        endIndex = 'endRowPlusOne'
        numIndices = 'numRows'
    else:
        raise ValueError ('Invalid sparse format "' + sparseFormat + '"')
    if dataLayout == 'row major':
        colRow = 'row'
        rowCol = 'column'
    elif dataLayout == 'column major':
        colRow = 'column'
        rowCol = 'row'
    else:
        raise ValueError ('Unknown data layout "' + dataLayout + '"')
    if conjugateMatrixEntries:
        briefConj = ',\n///   using conjugate of sparse matrix elements'
    else:
        briefConj = ''

    substDict = {'sparseFormat': sparseFormat, 
                 'colRow': colRow, 'rowCol': rowCol,
                 'briefConj': briefConj, 'numIndices': numIndices,
                 'startIndex': startIndex, 'endIndex': endIndex,
                 'fmtColRow': fmtColRow, 'fmtRowCol': fmtRowCol}

    brief = '/// ${sparseFormat} sparse matrix-(multi)vector multiply\n'
    brief = brief + '///   with ${colRow}-major input / output (multi)vectors\n'
    if conjugateMatrixEntries:
        brief = brief + '///   using conjugate of sparse matrix entries\n'

    body = '''///
/// \\tparam Ordinal The type of indices used to access the entries of
///   the sparse and dense matrices.  Any signed or unsigned integer
///   type which can be used in pointer arithmetic with raw arrays 
///   will do.
/// \\tparam MatrixScalar The type of entries in the sparse matrix.
///   This may differ from the type of entries in the input/output
///   matrices.
/// \\tparam DomainScalar The type of entries in the input multivector Y.
///   This may differ from the type of entries in the output multivector X.
/// \\tparam RangeScalar The type of entries in the output multivector X.
///
/// \param numRows [in] Number of rows in the sparse matrix.
/// \param numCols [in] Number of columns in the sparse matrix.
'''
    if not defDict['hardCodeNumVecs']:
        body = body + '''/// \param numVecs [in] Number of columns in X or Y (must be the same
///   for both).
'''
    body = body + '''/// \param X [out] Output multivector, stored in ${colRow}-major order.
/// \param LDX [in] Stride between ${colRow}s of X.  We assume unit
///   stride between ${rowCol}s of X.
/// \param ptr [in] Length (${numIndices}+1) array of index offsets 
///   between ${fmtRowCol}s of the sparse matrix.
/// \param ind [in] Array of ${fmtColRow} indices of the sparse matrix.
///   ind[ptr[i] .. ptr[i+1]-1] are the ${fmtColRow} indices of
///   ${fmtRowCol} i (zero-based) of the sparse matrix.
/// \param val [in] Array of entries of the sparse matrix.
///   val[ptr[i] .. ptr[i+1]-1] are the entries of ${fmtRowCol} i
///   (zero-based) of the sparse matrix.
/// \param Y [in] Input multivector, stored in ${colRow}-major order.
/// \param LDY [in] Stride between ${colRow}s of Y.  We assume unit
///   stride between ${rowCol}s of Y.'''

    doc = Template(brief + body + '\n').substitute (substDict)
    brief = '' # Release memory we don't need anymore    
    body = '' # Release memory we don't need anymore    
    # Indent each line.
    return '\n'.join(' '*indent + line for line in doc.split('\n'))

def emitHeaderDeclFile (filename):
    '''Make a header file with declarations of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of declarations for the sparse triangular solve
    routines.
    '''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice ()
    s = s + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Declarations of "raw" sequential sparse matrix-vector multiply routines.
/// \warning This code was generated by the SparseMatVec.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s = s + emitFuncDeclVariants(0) + \
        '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __' + headerizedFilename + '\n'
    return s

def emitHeaderDefFile (filename):
    '''Emit a header file with definitions of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of definitions for the sparse triangular solve
    routines.
    '''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice () + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Definitions of "raw" sequential sparse triangular solve routines.
/// \warning This code was generated by the SparseTriSolve.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s = s + emitFuncDefVariants(0) + \
        '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __' + headerizedFilename + '\n'
    return s

def run ():
    '''Generate the two header files mentioned in the module's documentation.

    This writes the header file of function declarations
    'Kokkos_Raw_SparseTriangularSolve_decl.hpp', and the header file
    of function definitions
    'Kokkos_Raw_SparseTriangularSolve_def.hpp', for all variants of
    sparse triangular solve that this module knows how to generate.
    Both files are written to the current working directory.
    '''

    rootName = 'Kokkos_Raw_SparseMatVec'
    declName = rootName + '_decl.hpp'
    defName = rootName + '_def.hpp'

    # No side effects (files opened for output or modified) until all
    # code generation has completed successfully.
    declStr = emitHeaderDeclFile (declName)
    defStr = emitHeaderDefFile (defName)

    # Write the header files.
    with open(declName, 'w') as declFile:
        declFile.write (declStr)
    with open(defName, 'w') as defFile:
        defFile.write (defStr)

# Code to execute if running the module as an executable script.
if __name__ == "__main__":
    import sys

    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else:
        run ()




