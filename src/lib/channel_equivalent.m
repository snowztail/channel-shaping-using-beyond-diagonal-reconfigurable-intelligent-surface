function [H_d, H_f, H_b] = channel_equivalent(H_d, H_f, H_b, W, G)
	H_d = pagemtimes(permute(G, [1, 2, 4, 3]), pagemtimes(H_d, W));
	H_f = pagemtimes(H_f, W);
	H_b = pagemtimes(G, H_b);
end
