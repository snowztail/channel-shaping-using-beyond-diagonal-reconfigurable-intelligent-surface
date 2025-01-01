function [W] = precoder_wsr(H, W, P_t, P_n, rho)
	[N_r, N_t, N_e, K] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4));
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-3, 0);
	iter.J = sum(rho .* rate_mimo(H, W, P_n), 3);
	while ~iter.converge && iter.counter <= 1e2
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		F = pagemtimes(H(:, :, logical(eye(K))), pageswap(W));
		G = pagemrdivide(pagectranspose(F), Q + T(:, :, logical(eye(K))));
		V = rho .* pagehermitize(eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));
		Y = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		for k = 1 : K
			fun = @(lambda) norm(mldivide(sum(Y(:, :, :, k), 3) + lambda * eye(N_t), H(:, :, k, k)' * G(:, :, k)' * V(:, :, k)), 'fro') ^ 2 - P_t;
			[lb, ub] = deal(0, 1e3);
			if fun(lb) * fun(ub) < 0
				lambda = bisection(fun, lb, ub);
			else
				lambda = 0;
			end
			W(:, :, :, k) = mldivide(sum(Y(:, :, :, k), 3) + lambda * eye(N_t), H(:, :, k, k)' * G(:, :, k)' * V(:, :, k));
		end
		J = sum(rho .* rate_mimo(H, W, P_n), 3);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end

function [A] = pagehermitize(A)
	A = 0.5 * (A + pagectranspose(A));
end

function [x] = bisection(fun, lb, ub)
	assert(lb < ub, 'Lower bound must be smaller than upper bound.');
	assert(fun(lb) * fun(ub) < 0, 'Initial values must have opposite signs.');

	epsilon = 1e-6;
    while (ub - lb) / 2 > epsilon
        x = (lb + ub) / 2;
        if fun(x) == 0
            break
        elseif fun(x) * fun(lb) < 0
            ub = x;
        else
            lb = x;
        end
    end
end
