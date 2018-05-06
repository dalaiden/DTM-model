function[conc_glace_interpol] = interpol_conc_glace_V3(x_model_traj, y_model_traj, data_for_interpol, type_interpo)

% Cette fonction donne la concentration de la glace a la position actuelle
% en fonction des donnees dans le voisinage (interpolation avec griddata)

% Calculer la concentration de glace pour chaque pas de temps a l'aide de la fonction griddata
conc_glace_interpol = griddata(data_for_interpol(:,4), data_for_interpol(:,5), data_for_interpol(:,8), ...
    x_model_traj, y_model_traj, type_interpo);

% Si la concentration en glace interpolee est superieur a 100 -> conc_glace_interpol = 100
conc_glace_interpol(conc_glace_interpol > 100) = 100;
 
% Si la concentration en glace interpolee est inferieur a 0 -> conc_glace_interpol = 0
conc_glace_interpol(conc_glace_interpol < 0) = 0;
end
