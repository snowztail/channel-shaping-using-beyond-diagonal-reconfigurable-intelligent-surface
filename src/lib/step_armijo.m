function [mu] = step_armijo(fun, Theta, D, S)
	O = fun(Theta);
	mu = 1;
	T = eye(size(Theta));
	T(S, S) = expm(mu * D(S, S));
	
	% * Undershoot, double the step size
	while (fun(T ^ 2 * Theta) - O) >= (mu * 0.5 * trace(D(S, S) * D(S, S)'))
		mu = mu * 2;
		T(S, S) = T(S, S) ^ 2;
	end
	
	% * Overshoot, halve the step size
	while (fun(T * Theta) - O) < (0.5 * mu * 0.5 * trace(D * D')) && (mu >= eps)
		mu = mu * 0.5;
		T(S, S) = expm(mu * D(S, S));
	end
end
