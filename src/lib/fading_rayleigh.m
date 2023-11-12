function [H] = fading_rayleigh(varargin)
	H = sqrt(0.5) * (randn(varargin{:}) + 1i * randn(varargin{:}));
end
