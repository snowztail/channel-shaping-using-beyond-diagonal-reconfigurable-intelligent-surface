function [B] = keep_block_diagonal(A, n)
% KEEP_BLOCK_DIAGONAL Keep block diagonal (size n x n) of square matrix A (size sn x sn) and set zero elsewhere
%   B = KEEP_BLOCK_DIAGONAL(A, n) returns a matrix B that is the same size as A, but with all elements outside the block diagonal of size n x n set to zero.

% Get size of A
[sn, ~] = size(A);

% Initialize B with zeros
B = zeros(sn);

% Set block diagonal of size n x n to A's block diagonal
for i = 1 : n : sn
    B(i : i + n - 1, i : i + n - 1) = A(i : i + n - 1, i : i + n - 1);
end

end
