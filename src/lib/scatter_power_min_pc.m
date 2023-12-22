function [Theta, H] = scatter_power_min_pc(H_d, H_f, H_b, L)
	persistent iter;
	if isempty(iter)
		Theta = eye(size(H_f, 1));
	else
		Theta = iter.Theta;
	end

	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[iter.converge, iter.tolerance, iter.counter, iter.P] = deal(false, 1e-4, 0, norm(H, 'fro') ^ 2);
	while ~iter.converge
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			D = H_d + H_b(:, S_c) * Theta(S_c, S_c) * H_f(S_c, :);
			T = H_b(:, S)' * H_b(:, S);
			B = max(eig(T)) * eye(L) - T + eps;
			M = (B * Theta(S, S) * H_f(S, :) - H_b(:, S)' * D) * H_f(S, :)';
			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		P = norm(H, 'fro') ^ 2;
		iter.converge = (abs(P - iter.P) / iter.P <= iter.tolerance);
		iter.P = P;
		iter.counter = iter.counter + 1;
	end
	iter.Theta = Theta;
end
