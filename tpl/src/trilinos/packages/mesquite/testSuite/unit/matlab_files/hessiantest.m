H = zeros(12,12)
b0 = [4 4 7; 4 5 7; 7 7 3]
b1 = [4 8 7; 3 5 7;1 2 3]
b2 = [4 4 7; 6 5 9; 1 8 5]
b3 = [4 4 2; 4 5 3; 2 3 3]
b4 = [2 4 7; 3 2 7; 1 4 3]
b5 = [8 4 9; 4 5 7; 9  7 3]
c0  = b0
c1 = b1
c2 = b2 + b5
c3 = b3
c4 = b4
c5 = b3
for i=1:3, for j=1:3, H(1+i-1, 1+j-1) = H(1+i-1, 1+j-1) + b0(i,j); end; end;
for i=1:3, for j=1:3, H(1+i-1, 4+j-1) = H(1+i-1, 4+j-1) + b1(i,j); end; end;
for i=1:3, for j=1:3, H(1+i-1, 7+j-1) = H(1+i-1, 7+j-1) + b2(i,j); end; end;
for i=1:3, for j=1:3, H(4+i-1, 4+j-1) = H(4+i-1, 4+j-1) + b3(i,j); end; end;
for i=1:3, for j=1:3, H(4+i-1, 7+j-1) = H(4+i-1, 7+j-1) + b4(i,j); end; end;
for i=1:3, for j=1:3, H(7+i-1, 7+j-1) = H(7+i-1, 7+j-1) + b5(i,j); end; end;

for i=1:3, for j=1:3, H(1+i-1, 1+j-1) = H(1+i-1, 1+j-1) + c0(i,j); end; end;
for i=1:3, for j=1:3, H(1+i-1, 10+j-1) = H(1+i-1, 10+j-1) + c1(i,j); end; end;
for i=1:3, for j=1:3, H(1+i-1, 4+j-1) = H(1+i-1, 4+j-1) + c2(i,j); end; end;
for i=1:3, for j=1:3, H(10+i-1, 10+j-1) = H(10+i-1, 10+j-1) + c3(i,j); end; end;
c4 = c4';
for i=1:3, for j=1:3, H(4+i-1, 10+j-1) = H(4+i-1, 10+j-1) + c4(i,j); end; end;
for i=1:3, for j=1:3, H(4+i-1, 4+j-1) = H(4+i-1, 4+j-1) + c5(i,j); end; end;
H
for i=1:12, for j=1:12, if (i>j) H(i,j) = H(j,i); end; end; end; 
H
x = [4 5 6 2 5 9 1 2 6 1 5 9]
H*x'
y = [3 2 6 1 2 4 3 6 9 2 4 4]
H*x'+y'



