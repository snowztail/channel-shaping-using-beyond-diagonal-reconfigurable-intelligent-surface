function [D] = direction_conjugate(G_r, iter)
	% * Periodic conjugate direction restart by steepest ascent
	if mod(iter.counter, numel(G_r)) == 0
		gamma = 0;
	else
		gamma = trace((G_r - iter.G_r) * G_r') / trace(iter.G_r * iter.G_r');
		% * If the conjugate direction is not ascent, restart
		if real(trace((G_r + gamma * iter.D)' * G_r)) < 0
			gamma = 0;
		end
	end
	D = G_r + gamma * iter.D;
end
