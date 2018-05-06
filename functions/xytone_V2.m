function [lat, lon] = xytone(x, y)

% Cette fonction peremt de convertir les coordonnees cartesiennes en
% coordonnees geographiques (projection des observations de vitesses de
% glace)

lon = nan(length(x),1);
lat = nan(length(x),1);

for i = 1:length(x);
    lon(i) = atan(y(i)/x(i));
    if sign(lon(i)) ~= sign(y(i));
        lon(i) = lon(i) + pi;
    end
    % Conversion en degree
    lon(i) = rad2deg(lon(i));
    lat(i) = 90 - (x(i)*0.001/(110.949*cos(deg2rad(lon(i)))));
    
    % Si longitude est inferieure a 0 degree -> + 360
    if lon(i) < 0;
        lon(i) = lon(i) + 360;
    end
end

end
