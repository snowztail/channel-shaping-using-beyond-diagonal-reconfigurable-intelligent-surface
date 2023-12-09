function [W] = precoder_rate_ic_bisection(H, W, P_t, P_n, rho)
	[N_r, N_t, N_e, K, P_t] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4), pagenorm(W, 'fro') .^ 2);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		F = pagemtimes(H(:, :, logical(eye(K))), pageswap(W));
		G = pagemrdivide(pagectranspose(F), Q + T(:, :, logical(eye(K))));
		V = shiftdim(rho, -2) .* (eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));
		Y = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));

		lambda = zeros(1, 1, 1, K);
		for k = 1 : K
			% fun = @(lambda_k) norm((sum(Y(:, :, :, k), 3) + lambda_k * eye(N_t)) \ H(:, :, k, k)' * G(:, :, k)' * V(:, :, k), 'fro') - P_t(k);
			fun = @(lambda_k) norm(inv((sum(Y(:, :, :, k), 3) + lambda_k * eye(N_t))) * H(:, :, k, k)' * G(:, :, k)' * V(:, :, k), 'fro') - P_t(k);
			if fun(0) * fun(1e12) > 0
				lambda(k) = 0;
			else
				lambda(k) = bisection(fun, 0, 1e12, 1e-6);
			end
			W(:, :, :, k) = inv(sum(Y(:, :, :, k), 3) + lambda(k) * eye(N_t)) * H(:, :, k, k)' * G(:, :, k)' * V(:, :, k);
		end

		Y = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		Z = pagemtimes(pagemtimes(W, 'ctranspose', Y, 'none'), W);
		lambda1 = real(shiftdim(pagetrace(P_n * pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G)) + pagetrace(sum(Z, 4)), -1) - pagetrace(sum(Z, 3))) ./ P_t;
		% % lambda = max(0, lambda);
		% X(:, :, 1, :) = pagemtimes(pagemtimes(G, H(:, :, logical(eye(K)))), 'ctranspose', V, 'none');
		% W = pagemldivide(sum(Y, 3) + lambda .* eye(N_t), X);
		J = rho' * rate_mimo_ic(H, W, P_n)
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		% iter.converge = all(abs(pagenorm(W, 'fro') .^ 2 - P_t) ./ P_t < iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
end

function [b] = pagetrace(A)
	b = sum(pageeig(A), 1);
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end

function [x] = bisection(fun, lb, ub, epsilon)
	assert(epsilon > 0, 'Tolerance must be non-negative.');
	assert(lb < ub, 'Lower bound must be smaller than upper bound.');
	assert(fun(lb) * fun(ub) < 0, 'Initial values must have opposite signs.');

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
    % x = (lb + ub) / 2;
end
