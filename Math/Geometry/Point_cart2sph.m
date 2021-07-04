function point = Point_cart2sph(point)

    point.r = sqrt(point.x.^2 + point.y.^2 + point.z.^2);
    
    point.theta = acos(point.z./point.r);
    point.theta(isnan(point.theta)) = 0;
    
    point.phi = atan2(point.y, point.x);
    
end
