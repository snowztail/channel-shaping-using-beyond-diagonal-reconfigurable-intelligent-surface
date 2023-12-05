function [H] = fading_los(N_r, N_t)
	H = steering_ula(N_r) * steering_ula(N_t)';
end

function [s] = steering_ula(N)
	azimuth = 2 * pi * rand;
	x = (0 : N - 1)';
	s = exp(-1i * pi * x * cos(azimuth));
end

function [s] = steering_upa(N_x, N_y)
	[azimuth, elevation] = deal(2 * pi * rand, 2 * pi * rand);
	[x, y] = meshgrid(0 : N_x - 1, 0 : N_y - 1);
	r = vec(x) * cos(azimuth) + vec(y) * sin(azimuth);
	s = exp(-1i * pi * r * cos(elevation));
end
