function [ B ] = normalizeMatrix( A )
%NORMALIZEMATRIX normalize the matrix
% 1) sum of absolute values equals 1
% 2) the spectral normalization
% 3) the value A(3,3) = 1

% 1) sum of the values equals 1
B = A * (1/sqrt(A(:)'*A(:)));

% % 2) the spectral normalization
% B = A * (1/(det(A)^(1/size(A,1))));
% 
% % 3) the value A(3,3) = 1
% B = A * (1/A(3,3));

end

