function [W] = precoder_leakage_ic(H, G)
    W = reciprocal(combiner_leakage_ic(reciprocal(H), reciprocal(G)));
end

function [Y] = reciprocal(X)
    Y = permute(pagectranspose(X), [1, 2, 4, 3]);
end
