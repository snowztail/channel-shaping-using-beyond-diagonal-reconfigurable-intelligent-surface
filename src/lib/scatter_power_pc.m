function [Theta, H] = scatter_power_pc(H_d, H_f, H_b, Theta, L)
	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	iter.P = norm(H, 'fro') ^ 2;
	while ~iter.converge
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			M = H_b(:, S)' * (H_d + H_b(:, S) * Theta(S, S) * H_f(S, :)) * H_f(S, :)';
			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		P = norm(H, 'fro') ^ 2;
		iter.converge = (abs(P - iter.P) / iter.P <= iter.tolerance);
		iter.P = P;
		iter.counter = iter.counter + 1;
	end
end
