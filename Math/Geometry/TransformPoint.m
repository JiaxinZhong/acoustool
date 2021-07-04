% 计算三维点在各种坐标系下的值

function point = TransformPoint(point)

	if ~isfield(point, 'x') 
		if isfield(point, 'r')
			point = Point_sph2cart(point);
		end
	end

	if ~isfield(point, 'r')
		if isfield(point, 'x')
			point = Point_cart2sph(point);
		end
	end

end
