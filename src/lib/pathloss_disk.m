function [L_d, L_f, L_b] = pathloss_disk(K, r, L, a_d, a_f, a_b)
	[C_t, C_s, C_r] = deal(distribution_disk(K, r)', 0, distribution_disk(K, r));
	[L_d, L_f, L_b] = deal(L * abs(C_r - C_t) .^ (-a_d), L * abs(C_s - C_t) .^ (-a_f), L * abs(C_r - C_s) .^ (-a_b));
end

function [C] = distribution_disk(K, r)
	C = sqrt(r ^ 2 * rand(K, 1)) .* exp(1i * 2 * pi * rand(K, 1));
end
