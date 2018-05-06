function [x, y] = netoxy(lon, lat)

% Cette fonction permt de convertir les coordonnees geographiques en
% coordonnees cartesiennes (projection des observations de vitesses de
% glace)

x = nan(length(lon),1);
y = nan(length(lon),1);

for i = 1:length(lon);
    x(i) = 110.949*(90 - lat(i))*cos(deg2rad(lon(i)));
    y(i) = 110.949*(90 - lat(i))*sin(deg2rad(lon(i)));
end

% Conversion en m et non en kilometres

x = x*1000;
y = y*1000;

end





