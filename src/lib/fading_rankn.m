function [H] = fading_rankn(N_r, N_t, N)
	H = zeros(N_r, N_t, N);
	for n = 1 : N
		H(:, :, n) = fading_los(N_r, N_t);
	end
	H = sqrt(1 / N) * sum(H, 3);
end
