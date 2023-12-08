function [W] = precoder_rate_ic(H, W, P_n, rho)
	[N_r, N_t, N_e, K, P_t] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4), shiftdim(pagenorm(W, 'fro')));

	W_ = W;

	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		F = pagemtimes(H(:, :, logical(eye(K))), W(:, :, 1 : end));
		G = pagemrdivide(pagectranspose(F), Q + T(:, :, logical(eye(K))));
		V = shiftdim(rho, -2) .* (eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));

		Z = pagemtimes(pagemtimes(pagemtimes(pagemtimes(G, H), W), 'ctranspose', V, 'none'), pagemtimes(pagemtimes(G, H), W));

		C1 = pagetrace(P_n * pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G));
		C2 = pagetrace(sum(Z, 4));
		C3 = pagetrace(sum(Z, 3));
		lambda_ = shiftdim(pagetrace(P_n * pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G)) + pagetrace(sum(Z, 4)), -1) - pagetrace(sum(Z, 3));

		% lambda_ = pagetrace(P_n * pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G))

		% O = pagemtimes(V, pagemtimes(pagemtimes(pagemtimes(G, H), W), 'none', pagemtimes(pagemtimes(G, H), W), 'ctranspose'));
		A = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		B = permute(pagemtimes(pagemtimes(V, G), F), [1, 2, 4, 3]);
		% lambda = 2 * pagetrace(pagemtimes(pagemtimes(W, 'ctranspose', sum(A, 3), 'none'), W)) - 2 * real(pagetrace(B));
		for k = 1 : K
			lambda(1, 1, 1, k) = trace(V(:, :, k) * G(:, :, k) * Q(:, :, k) * G(:, :, k)' - W(:, :, k)' * (sum(A(:, :, 1 : end ~= k, k), 3) * W(:, :, k)));
		end
		W = pagemldivide(sum(A, 3) + lambda .* eye(N_t), permute(pagemtimes(pagemtimes(G, H(:, :, logical(eye(K)))), 'ctranspose', V, 'none'), [1, 2, 4, 3]));

		norm(W(:, :, :, 1), 'fro') ^ 2
		J = rho' * rate_mimo_ic(H, W, P_n);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;


		% lambda = sum(pageeig(pagemtimes(pagemtimes(pagemtimes(V, G), Q), pagectranspose(G)) - (permute(sum(O, 3), [1, 2, 4, 3]) - O(:, :, logical(eye(K))))), 1);
		% lambda = sum(real(pageeig(pagemtimes(pagemtimes(pagemtimes(V, G), Q), pagectranspose(G)) - (permute(sum(O, 3), [1, 2, 4, 3]) - O(:, :, logical(eye(K)))))), 1);
		% W =
		% A = pagemtimes(pagemtimes(G, H), 'none', pagemtimes(G, H), 'ctranspose');
		% A = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		% W = pagemldivide(permute(sum(A, 3), [1, 2, 4, 3]) + lambda .* eye(N_t), pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'));
		% A = pagemtimes(pagemtimes(G, H), 'none', pagemtimes(G, H), 'ctranspose');
	end










	W = W_;
	P_t = norm(W(:, :, :, 1), 'fro') ^ 2;
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
	G = zeros(N_e, N_r, K);
	V = zeros(N_e, N_e, K);
	lambda = zeros(K, 1);
	lambda_ = zeros(K, 1);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		for k = 1 : K
			G(:, :, k) = W(:, :, :, k)' * H(:, :, k, k)' / (Q(:, :, k) + T(:, :, k, k));
			V(:, :, k) = rho(k) * (eye(N_e) + W(:, :, :, k)' * H(:, :, k, k)' / Q(:, :, k) * H(:, :, k, k) * W(:, :, :, k));
			% V(:, :, k) = (W(:, :, :, k)' * H(:, :, k, k)' / Q(:, :, k) * H(:, :, k, k) * W(:, :, :, k));
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
			% p = max(lambda(k), 0);
			p = lambda_(k);
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

function [b] = pagetrace(A)
	b = sum(pageeig(A), 1);
end
