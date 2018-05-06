function[rand_pts_canada] = formate_random_pts_nb(echantillon, cote_ORCA, distance_limite,nb_pts)

% Charger les donnees
data = importdata(echantillon);
data = data.data;

% La limite de la distance (zone ou les arbres partiront) -> distance_limite

% Transformer les lon : 0 -> 360
for i = 1:length(data);
    if data(i,2) < 0;
        data(i,2) = data(i,2) + 360;
    end
    % Ajouter la position en x/y (systeme de coordonnees des bouees)
    [data(i,4), data(i,5)] = netoxy(data(i,2), data(i,3));
    
    % Ajouter la distance par rapport a la cote la plus proche
    data(i,6) = distance_cote(cote_ORCA, data(i,4:5));
end

data = data(data(:,3) > 70,:);
data = data((data(:,6) > distance_limite(1)*10^3) & (data(:,6) < distance_limite(2)*10^3),:);

% Sorties
% 1/ LON
% 2/ LAT

rand_pts_canada = data(:,2:3);

rand_pts_canada = datasample(rand_pts_canada,nb_pts);

end
