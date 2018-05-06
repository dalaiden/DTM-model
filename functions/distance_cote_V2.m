function[coast_distance_vecteur] = distance_cote_V2(cote_ORCA, x_model_traj, y_model_traj, nb_arbres)

% Cette fonction donne comme output la distance du point actuel par rapport
% a la cote (V2 car pour tous arbres)
% ATTENTION, nb_arbres EST LE NOMBRE D'ARBRES QUE L'ON VEUT CALCULER. DONC
% SI ON TOURNE SUR UN ARBRE, nb_arbres = 1 !

% Premierement, il faut transformer les coordonnes des cotes en x/y
[x_y_cote(:,1), x_y_cote(:,2)] = netoxy(cote_ORCA(:,1), cote_ORCA(:,2));


coast_distance_vecteur = nan(1,nb_arbres);
for k = 1:nb_arbres

	distance_cote_matrix = nan(length(cote_ORCA),1);
	for w = 1:length(cote_ORCA);
		distance_x = x_model_traj(k) - x_y_cote(w,1);
		distance_y = y_model_traj(k) - x_y_cote(w,2);
		distance_cote_matrix(w) = sqrt((distance_x^2) + (distance_y^2));
	end
	% Pour savoir, si on est pres de la cote, prendre le minimum de la matrice distance cote
	coast_distance = min(distance_cote_matrix);

    % Stocker la distance minimale pour chaque arbre
    coast_distance_vecteur(k) = coast_distance;
 
    % Clean up 
    clear coast_distance
end
