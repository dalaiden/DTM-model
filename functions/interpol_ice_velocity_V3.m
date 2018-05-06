function[n,e] = interpol_ice_velocity_V3(x_model_traj, y_model_traj, data_for_interpol, type_interpo)

% Cette fonction interpole la vitesse de glace pour la position actuelle

n = griddata(data_for_interpol(:,4), data_for_interpol(:,5), data_for_interpol(:,7), ...
      x_model_traj, y_model_traj, type_interpo);
e = griddata(data_for_interpol(:,4), data_for_interpol(:,5), data_for_interpol(:,6), ...
      x_model_traj, y_model_traj, type_interpo);
end


