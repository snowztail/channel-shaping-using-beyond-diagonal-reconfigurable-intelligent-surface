function [S] = simplex_standard(K, delta)
	N = ceil(1 / delta);
	T = nchoosek(1 : K + N - 1, N) - (0 : N - 1);
	S = histc(T', 1 : K) / N;
end
