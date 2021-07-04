function point = cal_cyl2sph(point)

	point = Point_cyl2cart(point);
	point = cal_cart2sph(point);
end
