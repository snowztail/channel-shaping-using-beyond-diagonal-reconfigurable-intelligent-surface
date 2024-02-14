function [H] = fading_ricean(N_r, N_t, K)
	if isinf(K)
		H = fading_los(N_r, N_t);
	else
		H = sqrt(K / (1 + K)) * fading_los(N_r, N_t) + sqrt(1 / (1 + K)) * fading_nlos(N_r, N_t);
	end
end
