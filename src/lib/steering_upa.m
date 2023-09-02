function [steering] = steering_upa(nAnts, angle)
	arguments
		nAnts;
		angle.azimuth = 2 * pi * rand;
		angle.elevation = 2 * pi * rand;
	end
	[coordinate.x, coordinate.y] = meshgrid(0 : sqrt(nAnts) - 1);
	coordinate.r = coordinate.x(:) * cos(angle.azimuth) + coordinate.y(:) * sin(angle.azimuth);
	steering = exp(-1i * pi * coordinate.r * cos(angle.elevation));
end
