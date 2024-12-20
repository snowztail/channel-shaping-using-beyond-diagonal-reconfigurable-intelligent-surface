function [R] = rate_mimo(H, W, P_n)
	if size(H, 3) == 1
		R = log(det(eye(size(W, 2)) + (H * W)' * (H * W) / P_n));
	else
		[N_r, K] = deal(size(H, 1), size(H, 3));
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		X = eye(N_r) + pagemldivide(Q, T(:, :, logical(eye(K))));
		R = log(pagedet(X));
	end
	R = real(R);
end

function [d] = pagedet(A)
	d = prod(pagesvd(A), 1);
end
