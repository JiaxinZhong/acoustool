function point = Point_cyl2cart(point)

	point.x = point.rho .* cos(point.phi);
	point.y = point.rho .* sin(point.phi);
