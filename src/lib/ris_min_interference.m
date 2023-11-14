function [Theta] = ris_min_interference(H_d, H_f, H_b, Theta, G)
	[K, L] = deal(size(H_d, 3), length(Theta) / G);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	I = interference_leakage(H)
	while ~iter.converge
		iter.I = I;
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			D = H_d + pagemtimes(permute(H_b(:, setdiff(1 : end, S), :), [1, 2, 4, 3]), pagemtimes(Theta(setdiff(1 : end, S), setdiff(1 : end, S)), H_f(setdiff(1 : end, S), :, :)));
			B = zeros(L, L, K);
			for k = 1 : K
				B(:, :, k) = 2 * max(eig(H_b(:, S, k)' * H_b(:, S, k))) * eye(L) - H_b(:, S, k)' * H_b(:, S, k) + eps;
			end
            M = 0;
			for j = 1 : K
				for k = setdiff(1 : K, j)
					M = M + B(:, :, k) * Theta(S, S) * H_f(S, :, j) * H_f(S, :, j)' - H_b(:, S, k)' * D(:, :, j, k) * H_f(S, :, j)';
				end
			end
			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		I = interference_leakage(H)
		iter.converge = (abs(I - iter.I) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end
