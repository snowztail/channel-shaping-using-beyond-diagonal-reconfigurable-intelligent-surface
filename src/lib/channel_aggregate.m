function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	switch ndims(H_d)
		case 2
			% Point-to-point channel
			H = H_d + H_b * Theta * H_f;
		case 4
			% Interference channel
			H = H_d + pagemtimes(pagemtimes(H_b, Theta), H_f);
		otherwise
			error('Invalid number of dimensions');
	end
end
