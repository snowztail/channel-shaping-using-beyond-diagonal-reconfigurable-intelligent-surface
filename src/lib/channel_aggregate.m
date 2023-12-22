function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	if size(H_d, 3) == 1
		H = channel_aggregate_pc(H_d, H_f, H_b, Theta);
	else
		H = channel_aggregate_ic(H_d, H_f, H_b, Theta);
	end
end

function [H] = channel_aggregate_pc(H_d, H_f, H_b, Theta)
	H = H_d + H_b * Theta * H_f;
end

function [H] = channel_aggregate_ic(H_d, H_f, H_b, Theta)
	H = H_d + pagemtimes(pagemtimes(H_b, Theta), H_f);
end
