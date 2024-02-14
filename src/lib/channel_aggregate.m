function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	H = H_d + H_b * Theta * H_f;
end
