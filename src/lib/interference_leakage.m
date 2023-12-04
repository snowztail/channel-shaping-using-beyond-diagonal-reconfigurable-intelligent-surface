function [I] = interference_leakage(H, W, G)
	K = size(W, 4);
	H = pagemtimes(pagemtimes(G, H), W);
	I = norm(H(:, :, ~logical(eye(K))), 'fro') ^ 2;
end
