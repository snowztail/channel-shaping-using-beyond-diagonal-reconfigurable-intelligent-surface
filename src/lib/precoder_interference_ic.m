function [W] = precoder_interference_ic(H, G, P_t)
	N_e = size(G, 1);
    W = sqrt(P_t / N_e) * reciprocal(combiner_interference_ic(reciprocal(H), reciprocal(G)));
end

function [Y] = reciprocal(X)
    Y = permute(pagectranspose(X), [1, 2, 4, 3]);
end
