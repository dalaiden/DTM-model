function[big_data_interpol_2] = formate_data_ice(vel_nord,vel_est,lon,lat,sea_ice_conc)

% This function gets data to interpolate trajectories

% Calculate cartesian coordonnees of the model
x_model_proj = nan(size(lon,1),size(lon,2));
y_model_proj = nan(size(lon,1),size(lon,2));
for i = 1:size(vel_nord,1);
    for j = 1:size(vel_nord,2);
        [x_model_proj(i,j), y_model_proj(i,j)] = netoxy(lon(i,j), lat(i,j));
    end
end

big_data_interpol = cell(size(vel_nord,3),1);
for i = 1:size(vel_nord,3);
    big_data_interpol{i}(:,1) = reshape(lon, size(lon,1)*size(lon,2),1);
    big_data_interpol{i}(:,2) = reshape(lat, size(lon,1)*size(lon,2),1);
    big_data_interpol{i}(:,3) = nan(size(lon,1)*size(lon,2),1);
    big_data_interpol{i}(:,4) = reshape(x_model_proj, size(lon,1)*size(lon,2),1); % x
    big_data_interpol{i}(:,5) = reshape(y_model_proj, size(lon,1)*size(lon,2),1); % y
    big_data_interpol{i}(:,6) = reshape(vel_est(:,:,i), size(lon,1)*size(lon,2),1); % EST
    big_data_interpol{i}(:,7) = reshape(vel_nord(:,:,i), size(lon,1)*size(lon,2),1); % NORD
    big_data_interpol{i}(:,8) = reshape(sea_ice_conc(:,:,i), size(lon,1)*size(lon,2),1);
end

% Select only data that located on the North Hemispher

big_data_interpol_2 = cell(size(vel_nord,3),1);
for i = 1:size(vel_nord,3);
    lignes_lat_sup_60 = find(big_data_interpol{i}(:,2) >= 60);
    big_data_interpol_2{i} = big_data_interpol{i}(lignes_lat_sup_60,:);
    % Clean up
    clear lignes_lat_sup_60
end

end
