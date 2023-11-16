function [W] = transmit_min_interference(H, G)
    W = reciprocal(receive_min_interference(reciprocal(H), reciprocal(G)));
end

function [Y] = reciprocal(X)
    Y = permute(pagectranspose(X), [1, 2, 4, 3]);
end
