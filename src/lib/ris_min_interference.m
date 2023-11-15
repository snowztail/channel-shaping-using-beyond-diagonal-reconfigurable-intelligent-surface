function [Theta] = ris_min_interference(H_d, H_f, H_b, Theta, G)
	[K, L] = deal(size(H_d, 3), length(Theta) / G);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	I = interference_leakage(H)
	while ~iter.converge
		iter.I = I;
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			D = H_d + pagemtimes(pagemtimes(H_b(:, S_c, :), Theta(S_c, S_c)), permute(H_f(S_c, :, :), [1, 2, 4, 3]));
			B = zeros(L, L, K);
			for k = 1 : K
				B(:, :, k) = 2 * max(eig(H_b(:, S, k)' * H_b(:, S, k))) * eye(L) - H_b(:, S, k)' * H_b(:, S, k);
			end
tic
            M = 0;
			for k = 1 : K
				for j = setdiff(1 : K, k)
					M = M + B(:, :, k) * Theta(S, S) * H_f(S, :, j) * H_f(S, :, j)' - H_b(:, S, k)' * D(:, :, k, j) * H_f(S, :, j)';
				end
			end
toc
tic
			Q = pagemtimes(pagemtimes(pagemtimes(B, Theta(S, S)), permute(H_f(S, :, :), [1, 2, 4, 3])) - pagemtimes(pagectranspose(H_b(:, S, :)), D), permute(pagectranspose(H_f(S, :, :)), [1, 2, 4, 3]));
			M_ = sum(Q(:, :, ~logical(eye(K))), [3, 4]);
toc
            % M_ = sum(Q, [3, 4]);

			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		I = interference_leakage(H)
		iter.converge = (abs(I - iter.I) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end
