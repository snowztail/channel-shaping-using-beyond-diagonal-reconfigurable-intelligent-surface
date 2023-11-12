function [Theta, H] = ris_min_interference(H_d, H_f, H_b, Theta, G, U, P_t, V, P_n)
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-10, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	I = interference_leakage(H_d, H_f, H_b, Theta);
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

function [I] = interference_leakage(H, U, V)

end

function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	K = size(H_d, 3);
	H = zeros(size(H_d));
	for k1 = 1 : K
		for k2 = 1 : K
			H(:, :, k1, k2) = H_d(:, :, k1, k2) + H_b(:, :, k1) * Theta * H_f(:, :, k2);
		end
	end
end
