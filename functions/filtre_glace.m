function[row_glace_inf_seuil, cell_data_arbres] = filtre_glace(cell_data_arbres,seuil_glace,nb_arbres)

% Cette fonction permet de trouver l'endroit ou la concentration en glace est inferieure au seuil

% Nombre de composantes du cell array
[n_compo, n_col] = size(cell_data_arbres);

row_glace_inf_seuil = nan(nb_arbres,1);
for i = 1:nb_arbres;
    % Lignes ou la concentration de glace est inferieure au seuil
    lines_seuil_inf = find((cell_data_arbres{5,2}(:,i) < seuil_glace) & (cell_data_arbres{2,2}(:,i) < 85));
    line_seuil_inf  = min(lines_seuil_inf);

    % On enregistre la ligne
    if isempty(line_seuil_inf) == 1;
       row_glace_inf_seuil(i) = nan;
    else
       row_glace_inf_seuil(i) = line_seuil_inf;
    end

    % On modifie toutes les composantes pour chaque arbre
    for j = 1:n_compo;
        [taille_mat, n_col_mat] = size(cell_data_arbres{j,2});
        if ((line_seuil_inf+1 > taille_mat) | (isempty(line_seuil_inf) == 1));
           cell_data_arbres{j,2}(:,i) = cell_data_arbres{j,2}(:,i);
        else  
           cell_data_arbres{j,2}(line_seuil_inf+1:taille_mat,i) = nan;
        end 
    end
end

end
