function [distance] = distance_2_pts(lon_1,lat_1,lon_2,lat_2)

[x_1, y_1] = netoxy(lon_1, lat_1);
[x_2, y_2] = netoxy(lon_2, lat_2);

distance_x = x_2 - x_1;
distance_y = y_2 - y_1;

distance = sqrt(distance_x^2 + distance_y^2);

% Convertion en km
distance = distance/1000;

end