function point = Point_sph2cart(point)

point.x = point.r .* cos(point.phi) .* sin(point.theta);
point.y = point.r .* sin(point.phi) .* sin(point.theta);
point.z = point.r .* cos(point.theta);
