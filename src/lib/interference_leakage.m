function [I] = interference_leakage(H)
	I = norm(H(:, :, ~logical(eye(size(H, 3)))), 'fro') ^ 2;
end
