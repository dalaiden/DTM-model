function[coast_distance] = distance_cote(cote_ORCA, x_y_model_traj)

% Cette fonction donne comme output la distance du point actuel par rapport
% a la cote

% Premierement, il faut transformer les coordonnes des cotes en x/y
[x_y_cote(:,1), x_y_cote(:,2)] = netoxy(cote_ORCA(:,1), cote_ORCA(:,2));

distance_cote_matrix = nan(length(cote_ORCA),1);
for w = 1:length(cote_ORCA);
    distance_x = x_y_model_traj(1) - x_y_cote(w,1);
    distance_y = x_y_model_traj(2) - x_y_cote(w,2);
    distance_cote_matrix(w) = sqrt((distance_x^2) + (distance_y^2));
end
% Pour savoir, si on est pres de la cote, prendre le minimum de la matrice distance cote
coast_distance = min(distance_cote_matrix);

end
