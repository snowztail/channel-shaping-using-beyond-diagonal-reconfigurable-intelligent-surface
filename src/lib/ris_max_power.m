function [Theta, H] = ris_max_power(H_d, H_f, H_b, Theta, G)
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-10, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	E = norm(H, 'fro') ^ 2;
	while ~iter.converge
		iter.E = E;
		for g = 1 : G
			S = (g - 1) / G * length(Theta) + 1 : g / G * length(Theta);
			M = H_b(:, S)' * (H_d + H_b(:, S) * Theta(S, S) * H_f(S, :)) * H_f(S, :)';
			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		E = norm(H, 'fro') ^ 2;
		iter.converge = (abs(E - iter.E) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end
