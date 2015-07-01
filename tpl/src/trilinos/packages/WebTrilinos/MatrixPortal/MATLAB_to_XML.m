%
% ierr = MATLAB_to_XML(FileName, A, LHS, RHS, ExactSolution)
%
% Writes a linear system in EpetraExt/XML format. The linear system is 
%   A * LHS = RHS
% A is a sparse matrix, and LHS and RHS are the left-hand side and the
% right-hand side.  ExactSolution is the exact solution, if known. All vectors
% can be multivectors (that is, have more than one column).
%
% The function returns a non-zero error code is an error occurs.
%
function ierr = MATLAB_to_XML(FileName, A, LHS, RHS, ExactSolution)

  n = size(A, 1);
  % problems should be square
  if n ~= size(A, 2)
    disp('the matrix is not square')
    ierr = -1
    return
  end
  if n ~= size(LHS, 1)
    disp('LHS is not compatible with A')
    ierr = -1
    return
  end
  if n ~= size(RHS, 1)
    disp('RHS is not compatible with A')
    ierr = -1
    return
  end
  m = size(LHS, 2);
  if m ~= size(RHS, 2)
    disp('LHS and RHS do not have the same number of vectors')
    ierr = -1
    return
  end
  
  fid = fopen(FileName, 'w');
  if fid == 0
    ierr = 03
    return
  end

  fprintf(fid, '<ObjectCollection Label="MATLAB generated">\n');

  % first write a Map
  fprintf(fid, '<Map Label="map" NumElements="%d" IndexBase="0" NumProc="1" ElementsOnProc0="%d">\n', n, n);
  fprintf(fid, '<Proc ID="0">\n');
  for i = 1:n
    fprintf(fid, '%d\n', i  - 1);
  end
  fprintf(fid, '</Proc>\n');
  fprintf(fid, '</Map>\n');

  % then the LHS vector

  fprintf(fid, '<MultiVector Label="LHS" Length="%d" NumVectors="%d" Type="double">\n', n, m);
  for i = 1:n
    for j =1:m
      fprintf(fid, '%e ', LHS(i, j));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '</MultiVector>\n');

  % then the RHS vector

  fprintf(fid, '<MultiVector Label="RHS" Length="%d" NumVectors="%d" Type="double">\n', n, m);
  for i = 1:n
    for j =1:m
      fprintf(fid, '%e ', RHS(i, j));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '</MultiVector>\n');

  % then the Exact vector

  fprintf(fid, '<MultiVector Label="ExactSolution" Length="%d" NumVectors="%d" Type="double">\n', n, m);
  for i = 1:n
    for j =1:m
      fprintf(fid, '%e ', ExactSolution(i, j));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '</MultiVector>\n');

  % then the matrix

  [rows,cols,vals] = find(A);
  nnz = size(rows, 1);
  
  fprintf(fid, '<PointMatrix Label="A" Rows="%d" Columns="%d" Nonzeros="%d" Type="double" StartingIndex="0">\n', n, n, nnz);
  for i = 1:nnz
    fprintf(fid, '%d %d %.15e\n', rows(i) - 1, cols(i) - 1, vals(i));
  end
  fprintf(fid, '</PointMatrix>\n');

  % finally close the file
  fprintf(fid, '</ObjectCollection>\n');

  fclose(fid);

  ierr = 0;
  return
