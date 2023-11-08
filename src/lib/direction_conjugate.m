function [D] = direction_conjugate(G_r, S, iter)
	L = length(S);
	
	% * Periodic conjugate direction restart by steepest ascent
	if mod(iter.counter, L ^ 2) == 0
		gamma = 0;
		iter.D(S, S) = zeros(L);
	else
		gamma = trace((G_r(S, S) - iter.G_r(S, S)) * G_r(S, S)') / trace(iter.G_r(S, S) * iter.G_r(S, S)');
		
		% * If the conjugate direction is not ascent, restart
		if real(trace((G_r(S, S) + gamma * iter.D(S, S))' * G_r(S, S))) < 0
			gamma = 0;
		end
	end
	D(S, S) = G_r(S, S) + gamma * iter.D(S, S);
end
