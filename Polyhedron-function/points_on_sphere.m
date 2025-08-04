function pts = points_on_sphere(N,ro)

gr = (sqrt(5.0) + 1.0) / 2.0;
ga=(2.0 - gr) * (2.0*pi);

pts = zeros(N,3);
for i = 1:N
    lat = asin(-1.0 + 2.0 * i / (N+1));
    lon = ga * i;
    tmp = [cos(lon)*cos(lat) sin(lon)*cos(lat) sin(lat)];
    tmp = ro*tmp;
    pts(i,:) = tmp;
end

end
