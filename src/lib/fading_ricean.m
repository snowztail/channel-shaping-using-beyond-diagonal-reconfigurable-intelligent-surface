function [H] = fading_ricean(N_t, type_t, N_r, type_r, K)
	arguments
		N_t;
		type_t;
		N_r;
		type_r;
		K = 0;
	end
	H = sqrt(K / (1 + K)) * component_los(steering(N_t, type_t), steering(N_r, type_r)) + sqrt(1 / (1 + K)) * component_nlos(N_t, N_r);
end

function [H] = component_los(s_t, s_r)
	H = s_r * s_t';
end

function [H] = component_nlos(N_t, N_r)
	H = sqrt(0.5) * (randn(N_r, N_t) + 1i * randn(N_r, N_t));
end

function [s] = steering(N, type)
	[azimuth, elevation] = deal(2 * pi * rand, 2 * pi * rand);
	
	switch type
		case 'ula'
			x = (0 : N - 1)';
			s = exp(-1i * pi * x * cos(azimuth));
		case 'upa'
			[N_x, N_y] = factor_pair_closest(N);
			[x, y] = meshgrid(0 : N_x - 1, 0 : N_y - 1);
			r = x(:) * cos(azimuth) + y(:) * sin(azimuth);
			s = exp(-1i * pi * r * cos(elevation));
	end
end

function [a, b] =  factor_pair_closest(n)
	amax = floor(sqrt(n));
	
	if mod(n, amax) == 0
		a = amax;
	else
		primeFactors  = factor(n);
		candidates = 1;
		for i = 1 : numel(primeFactors)
			f = primeFactors(i);
			candidates  = union(candidates, f .* candidates);
			candidates(candidates > amax) = [];
		end
		a = candidates(end);
	end
	
	b = n / a;
end
