function B = vec(A);

[n1,n2]=size(A);
B = reshape(A,n1*n2,1);
