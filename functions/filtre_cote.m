function[data_filtre_cote] = filtre_cote(cell_data_arbres,seuil_cote,nb_arbres)

% Cette fonction permet de ne garder que les positions des arbres dont la distance
% par rapport a la cote est inferieure au seuil.

% Nombre de composantes du cell array
[n_compo, n_col] = size(cell_data_arbres);
for i = 1:nb_arbres;
    % Lignes ou le seuil est superieure a la distance entre le point et la cote la plus proche
    lines_seuil_inf = find(cell_data_arbres{6,2}(:,i) < seuil_cote);
    line_seuil_inf  = min(lines_seuil_inf);

    % On modifie toutes les composantes pour chaque arbre
    for j = 1:n_compo;
        [taille_mat, n_col_mat] = size(cell_data_arbres{j,2});
        if taille_mat > line_seuil_inf;
           cell_data_arbres{j,2}(line_seuil_inf+1:taille_mat,i) = nan;
        end
    end
end

% Output
data_filtre_cote = cell_data_arbres;

end
