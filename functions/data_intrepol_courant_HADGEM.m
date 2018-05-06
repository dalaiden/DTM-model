function[data_for_interpol] = data_intrepol_courant_HADGEM(big_data_interpol, num_data)

% Cette fonction donne les observations des 8 points les plus proches
% georgraphiquement de la position actuelle.

% 2 inputs : 1/ Base de donnees mensuelles
%            2/ numero du jour (num_data)

% Etant donne que les donnees des courants oceaniques ne sont disponibles
% que tous les mois, il faut donc prendre la donnee la plus proche dans
% le temps. La variable temporelle (jour de l'annee) pour les donnees des 
% vitesses de courant se trouve a la deuxieme colonne du cell array. Pour
% trouver l'observation la plus proche, il faut donc minimiser le 
% abs(num_data - jour)

for i = 1:12;
    diff_jour(i,1) = abs(big_data_interpol{i,2} - num_data);
end
% Trouver le minimum
min_diff  = min(diff_jour);
num_ligne = find(min_diff == diff_jour);
num_ligne = num_ligne(1);

% Output
data_for_interpol = big_data_interpol{num_ligne,1};

end
