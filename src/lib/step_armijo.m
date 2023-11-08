function [Theta] = step_armijo(fun, Theta, D)
	O = fun(Theta);
	mu = 1;
	T = expm(mu * D);
	% * Undershoot, double the step size
	while (fun(T ^ 2 * Theta) - O) >= (mu * 0.5 * trace(D * D'))
		mu = mu * 2;
		T = T ^ 2;
	end
	% * Overshoot, halve the step size
	while (fun(T * Theta) - O) < (0.5 * mu * 0.5 * trace(D * D')) && (mu >= eps)
		mu = mu * 0.5;
		T = expm(mu * D);
	end
	Theta = T * Theta;
end
