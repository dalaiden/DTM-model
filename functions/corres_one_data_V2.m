function[lot_year_one_data] = corres_one_data(nb_annee,nb_days_data)

% Correspondance entre la ligne de l'iteration et la ligne
% de la base de donnee. Ici le probleme est qu'on fait tourner
% le modele sur 10 ans mais avec des donnees de la meme annee.
% Il faut donc que le programme comprenne qu'apres 365 (et ses
% multiples), il reprenne la premiere donnee disponibles.

% Input : nombre d'annees
% Outputs : 1/Iterasion
%           2/numero de donnee a reprendre 

nb_annee = nb_annee;
lot_year_one_data = nan(nb_annee*nb_days_data,2);

debut_fin_annee = nan(nb_annee,2);
debut_fin_annee(1,1) = 1;
debut_fin_annee(1,2) = nb_days_data;

for i = 2:nb_annee;
    debut_fin_annee(i,1) = debut_fin_annee(i-1,2) + 1;
    debut_fin_annee(i,2) = debut_fin_annee(i-1,2) + nb_days_data;
end

for i = 1:length(lot_year_one_data);
    lot_year_one_data(i,1) = i;
end

for i = 1:nb_annee;
    lot_year_one_data(debut_fin_annee(i,1):debut_fin_annee(i,2),2) = 1:nb_days_data;
end

end
