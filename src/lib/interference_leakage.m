function [I] = interference_leakage(H)
	K = size(H, 3);
	I = norm(H(:, :, ~logical(eye(K))), 'fro') ^ 2;
end
