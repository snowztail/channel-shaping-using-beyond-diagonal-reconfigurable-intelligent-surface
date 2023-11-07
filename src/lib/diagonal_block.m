function [B] = diagonal_block(A, n)
	B = zeros(size(A));
	for i = 1 : n : length(B)
		B(i : i + n - 1, i : i + n - 1) = A(i : i + n - 1, i : i + n - 1);
	end
end
