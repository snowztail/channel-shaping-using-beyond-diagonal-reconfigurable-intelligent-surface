function [W] = precoder_rate_ic(H, W, P_n, rho)
	[d, N_r, N_t, K] = deal(size(W, 2), size(H, 1), size(H, 2), size(H, 3));
	P_t = norm(W(:, :, :, 1), 'fro') ^ 2;
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
	G = zeros(d, N_r, K);
	V = zeros(d, d, K);
	lambda = zeros(K, 1);
	lambda_ = zeros(K, 1);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		for k = 1 : K
			G(:, :, k) = W(:, :, :, k)' * H(:, :, k, k)' / (Q(:, :, k) + T(:, :, k, k));
		end
		for k = 1 : K
			V(:, :, k) = rho(k) \ (eye(d) + W(:, :, :, k)' * H(:, :, k, k)' * Q(:, :, k) * H(:, :, k, k) * W(:, :, :, k));
		end
		for k = 1 : K
			A = 0;
			for j = 1 : K
				A = A + H(:, :, j, k)' * G(:, :, j)' * V(:, :, j) * G(:, :, j) * H(:, :, j, k);
			end
			lambda(k) = real(trace(V(:, :, k) * G(:, :, k) * Q(:, :, k) * G(:, :, k)') - trace(W(:, :, k)' * (A - H(:, :, k, k)' * G(:, :, k)' * V(:, :, k) * G(:, :, k) * H(:, :, k, k)) * W(:, :, k))) / P_t;
			B = -trace(V(:, :, k) * G(:, :, k) * P_n * eye(N_r) * G(:, :, k)');
			for l = setdiff(1 : K, k)
				B = B + trace(V(:, :, l) * G(:, :, l) * H(:, :, l, k) * W(:, :, k) * W(:, :, k)' * H(:, :, l, k)' * G(:, :, l)') ...
					  - trace(V(:, :, k) * G(:, :, k) * H(:, :, k, l) * W(:, :, l) * W(:, :, l)' * H(:, :, k, l)' * G(:, :, k)');
			end
			lambda_(k) = real(- B / P_t);
			p = max(lambda(k), 0);
			W(:, :, :, k) = (A + p * eye(N_t)) \ H(:, :, k, k)' * G(:, :, k)' * V(:, :, k);
			T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
			Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
			norm(W(:, :, :, k), 'fro') ^ 2
		end
		J = rho' * rate_mimo_ic(H, W, P_n);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
end
