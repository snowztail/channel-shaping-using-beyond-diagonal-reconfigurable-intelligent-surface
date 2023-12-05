function [Theta, H] = reflector_power_pc(H_d, H_f, H_b, Theta, L)
	G = length(Theta) / L;
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-10, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	P = norm(H, 'fro') ^ 2;
	while ~iter.converge
		iter.P = P;
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			M = H_b(:, S)' * (H_d + H_b(:, S) * Theta(S, S) * H_f(S, :)) * H_f(S, :)';
			[U, ~, V] = svd(M);
			Theta(S, S) = U * V';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		P = norm(H, 'fro') ^ 2;
		iter.converge = (abs(P - iter.P) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end
