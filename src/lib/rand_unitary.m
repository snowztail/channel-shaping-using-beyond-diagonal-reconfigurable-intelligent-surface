function [X] = rand_unitary(varargin)
	X = randn(varargin{:}) + 1i*randn(varargin{:});
	X = X * (X' * X) ^ -0.5;
end
