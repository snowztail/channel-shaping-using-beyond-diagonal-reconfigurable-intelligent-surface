function [W] = transmit_min_interference(H, G)
	% W = pagectranspose(receive_min_interference(permute(pagectranspose(H), [1, 2, 4, 3]), permute(pagectranspose(G), [1, 2, 4, 3])));
    % W = pagectranspose(receive_min_interference(permute(pagectranspose(H), [1, 2, 4, 3]), permute(pagectranspose(G), [1, 2, 4, 3])));
    W = reciprocal(receive_min_interference(reciprocal(H), reciprocal(G)));
end

function [Y] = reciprocal(X)
    Y = permute(pagectranspose(X), [1, 2, 4, 3]);
end
