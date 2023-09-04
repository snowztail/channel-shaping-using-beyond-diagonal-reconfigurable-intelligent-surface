function [R] = rate_mimo(H, Q, P_n)
	R = log(det(eye(size(H, 1)) + (H * Q * H') / P_n));
	R = real(R);
end
