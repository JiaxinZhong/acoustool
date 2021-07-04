function dist = DistanceBetweenPoints(point1, point2)

% 	point1 = transform_point(point1);
% 	point2 = transform_point(point2);
	dist = sqrt((point1.x-point2.x).^2 ...
		+ (point1.y-point2.y).^2 ...
		+ (point1.z-point2.z).^2);
end
