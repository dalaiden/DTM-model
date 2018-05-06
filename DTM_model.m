% DTM MODEL (Driftwood Transport Model)
%
% This model simulates Arctic driftwood trajectories
%
% 3 inputs :
%            1) Sea ice velocities    (lat,lon,time) DAILY
%            2) Sea ice concentration (lat,lon,time) DAILY
%            3) Sea surface currents  (lat,lon,time) MONTHLY
%
% lon and lat must be the same for the three inputs !
%
% Author : Quentin Dalaiden (quentin.dalaiden@uclouvain.be)
%
% Institute : Universite catholique de Louvain

function DTM_model(ice_vel_nord,...       % Sea ice velocity North-component
                   ice_vel_est,...        % Sea ice velocity East-component
                   LAT_M,...              % Latitude
                   LON_M,...              % Longitude
                   sic_M,...              % Sea ice concentrations
                   ocean_vel_nord,...     % Sea surface current North-component
                   ocean_vel_est,...      % Sea surface current East-component
                   duree_simul,...        % Time duration in years 
                   seuil_glace,...        % Ice threshold (0 to 1) 
                   seuil_cote,...         % Coast threshold in meters       
                   duree_derive_conc,...  % Time duration of the drift with sea surface current 
                   nb_arbres,...          % Number of woods
                   nb_waves,...           % Number of waves
                   name_EXP,...           % Name of the experiment
                   SHOW_FIG,...           % Show figures ?
                   PLOT_TRAJ_O,...        % Plot of each trajectory ?
                   day_ocean_data)        % day of the data ocean

%-----------------------------
%
%  BEGINNING OF THE CODE
%
%-----------------------------

warning off
addpath(genpath('functions/'))

%----------------------------
% PARAMETERS
%----------------------------

start_days         = sort(randi([90 270],1,nb_waves));
zone_lancement     = [300 400];
distance_exclusion = 1000;

%----------------------------

% Load coasts mask
load('data/cote_ORCA.mat')

%----------------------------
% Inputs formating
%----------------------------
disp('Working on loading data ...')

big_data_interpol     = formate_data_ice(ice_vel_nord,ice_vel_est,LON_M,LAT_M,sic_M); % ICE
data_courant_interpol = formate_data_ocean(ocean_vel_nord,ocean_vel_est,LON_M,LAT_M,day_ocean_data); % OCEAN 

for s = 1:length(start_days);

    % Clean up
    %clearvars -except LAT_M LON_M duree_simul seuil_glace seuil_cote duree_derive_conc nb_arbres nb_waves start_days zone_lancement distance_exclusion big_data_interpol data_courant_interpol cote_ORCA ice_vel_nord ice_vel_est sic_M ocean_vel_nord ocean_vel_est day_ocean_data s name_EXP PLOT_TRAJ_O SHOW_FIG

    clearvars -except LAT_M LON_M duree_simul seuil_glace seuil_cote duree_derive_conc nb_arbres nb_waves start_days zone_lancement distance_exclusion big_data_interpol data_courant_interpol cote_ORCA s name_EXP PLOT_TRAJ_O SHOW_FIG

    jour_debut = start_days(s);

    % Ici, il faut mettre les donnees de depart pour les arbres pour chaque region
    data_loc_canada  = formate_random_pts_nb('data/pts_canada_600_kil.csv', cote_ORCA, zone_lancement,nb_arbres);
    data_loc_siberie = formate_random_pts_nb('data/pts_siberie_600_kil.csv', cote_ORCA, zone_lancement,nb_arbres); 
    data_loc_arbres  = [data_loc_canada;data_loc_siberie];

    % Constuire la matrice qui contient le resume de la simulation
    data_recap = nan(nb_arbres,15);

    % Creer les dossiers de sortie
    namefolder = sprintf('outputs');
    mkdir(sprintf('%s/figures',namefolder))

    disp('Working on trajectory interpolation (ice velocity) ...')

    % Constuire la matrice des coordonnees des arbres
    cell_data_arbres = cell(1,1);
    % 1/ LON
    % 2/ LAT
    % 3/ X
    % 4/ Y
    % 5/ Concentraion glace
    cell_data_arbres{1,1} = 'LON';
    cell_data_arbres{2,1} = 'LAT';
    cell_data_arbres{3,1} = 'X';
    cell_data_arbres{4,1} = 'Y';
    cell_data_arbres{5,1} = 'Conc_glace';
    cell_data_arbres{6,1} = 'Distance_cote_nearest';

    % Mettre les coordonnes initiales   
    cell_data_arbres{1,2}(1,:) = transpose(data_loc_arbres(:,1));
    cell_data_arbres{2,2}(1,:) = transpose(data_loc_arbres(:,2));
    [cell_data_arbres{3,2}(1,:), cell_data_arbres{4,2}(1,:)] = netoxy(cell_data_arbres{1,2}(1,:), cell_data_arbres{2,2}(1,:));

    % Trajectory interpolation
    % When iteration (i) is higher than time size of data, take back the first data
    if size(big_data_interpol,1) > 360;
       correspondance_iterations_data = corres_data_nb_years();
       size_t_E = (5*size(big_data_interpol,1))-1;
       disp('INTERANNUAL')
    else  
      correspondance_iterations_data = corres_one_data_V2(duree_simul,size(big_data_interpol,1));
      size_t_E = (duree_simul*size(big_data_interpol,1))-1;
      disp('ANNUAL')
    end
    
    z = 1;
    % Faire tourner le modele pendant la duree de la simulation

    for i = jour_debut:size_t_E;
        % When iteration (i) is higher than time size of data, take back the first data
        num_data = correspondance_iterations_data(i,2);
        
        % Definir les coordonnees geographiques a chaque instant        
        [cell_data_arbres{2,2}(z,:), cell_data_arbres{1,2}(z,:)] = xytone_V2(cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:));

        % Definir les coordonnees cartesiennes a chaque instant
        [cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:)] = netoxy(cell_data_arbres{1,2}(z,:), cell_data_arbres{2,2}(z,:));

        % Donnees pour l'interpo
        data_for_interpol = big_data_interpol{num_data};
 
        % Calculer la vitesse de la glace
        [n,e] = interpol_ice_velocity_V3(cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:), data_for_interpol, 'linear');

        % Calculer la concentration de la glace
        cell_data_arbres{5,2}(z,:) = interpol_conc_glace_V3(cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:),...
                                     data_for_interpol, 'linear');

        % Calculer la distance par rapport a la cote la plus proche
        cell_data_arbres{6,2}(z,:) = distance_cote_V2(cote_ORCA, cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:), nb_arbres*2);

        % Convertir n et e en u et v
    	  v_model = e.*cos(deg2rad(cell_data_arbres{1,2}(z,:))) - (n.*sin(deg2rad(cell_data_arbres{1,2}(z,:))));
    	  u_model = -(n.*cos(deg2rad(cell_data_arbres{1,2}(z,:)))) - (e.*sin(deg2rad(cell_data_arbres{1,2}(z,:))));

    	  % Euler forward : interpolation de position
        cell_data_arbres{3,2}(z+1,:) = cell_data_arbres{3,2}(z,:) + (u_model.*86400);
        cell_data_arbres{4,2}(z+1,:) = cell_data_arbres{4,2}(z,:) + (v_model.*86400); 

        % Iterations
        z = z+1;
    end
    % Ajouter la conversion des x y en lon lat + calcul de la concenation de glace + distance cote pour la derniere ligne
    [cell_data_arbres{2,2}(z,:), cell_data_arbres{1,2}(z,:)] = xytone_V2(cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:));
    cell_data_arbres{5,2}(z,:) = interpol_conc_glace_V3(cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:),data_for_interpol, 'linear');
    cell_data_arbres{6,2}(z,:) = distance_cote_V2(cote_ORCA, cell_data_arbres{3,2}(z,:), cell_data_arbres{4,2}(z,:), nb_arbres*2);
 
    % Supprimer les jours ou l'arbre est trop pres de la cote (seuil cote);
    cell_data_arbres = filtre_cote(cell_data_arbres,seuil_cote,nb_arbres*2);

    % Localiser les points ou la concentration de glace est inferieure au seuil (pour chaque arbre).
    % A partir de cette ligne, supprimer le calcul de l'interpolation de trajectoires
    [row_glace_inf_seuil,cell_data_arbres] = filtre_glace(cell_data_arbres, seuil_glace,nb_arbres*2);

    % Tourner le modele avec les vitesses du courant oceanique pendant "duree_derive_conc"
    % Afin d'optimiser cette partie, ne prendre que les arbres qui ont au moins drifte pendant 365 jours
    % De plus, ne pas prendre les NaN
    ID_arbres_drift_ocean = 1:nb_arbres*2;
    ID_arbres_drift_ocean = transpose(ID_arbres_drift_ocean);
    ID_arbres_drift_ocean(:,2) = transpose(row_glace_inf_seuil); 
    % Connaitre les arbres a drifter avec les courants marins
    arbres_a_drifter = find((ID_arbres_drift_ocean(:,2) > 365) & (isnan(ID_arbres_drift_ocean(:,2)) == 0));
    ID_arbres_drift_ocean = ID_arbres_drift_ocean(arbres_a_drifter,:);
    %save('ID_arbres_drift_ocean.mat','ID_arbres_drift_ocean')
    disp('Working on trajectory interpolation (ocean''s courant velocity) ...')
    
    % Interpolation
    for i = 1:length(arbres_a_drifter);
        % Recommencer l'interpolation a partir du premier jour ou la concentration en glace est inferieure au seuil
        z = ID_arbres_drift_ocean(i,2);
        for t = 1:duree_derive_conc;
            % Connaitre le jour dans les donnees de la simulation
            jour_simu = z + (jour_debut-1);

            % When iteration (i) is higher than time size of data, take back the first data
            num_data = correspondance_iterations_data(jour_simu,2);

            if num_data <= size(big_data_interpol,1)*duree_simul;
		        % Definir les coordonnees geographiques a chaque instant
		        [cell_data_arbres{2,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1))] = xytone_V2(cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)));

		        % Matrice pour l'interpolation des vitesses des courants oceaniques. Vu que les donnees des courants oceaniques ne sont disponibles que
		        % tous les 5 jours, prendre la donnee la plus proche dans le temps         
		        data_for_interpol_courant = data_intrepol_courant_HADGEM(data_courant_interpol, num_data);

		        % Donnees pour l'interpo (concentration de la glace)
		        data_for_interpol = big_data_interpol{num_data};
		          
		        % Calculer la concentration de la glace
		        cell_data_arbres{5,2}(z,ID_arbres_drift_ocean(i,1)) = interpol_conc_glace_V3(cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)),...
		                                 data_for_interpol, 'linear');

		        % Calculer la distance par rapport a la cote la plus proche
		        cell_data_arbres{6,2}(z,ID_arbres_drift_ocean(i,1)) = distance_cote_V2(cote_ORCA, cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)), 1);

		        % Calculer la vitesse du courant oceanique
		        [n_courant,e_courant] = interpol_courant_velocity_V3(cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)), data_for_interpol_courant, 'linear');

		        % Convertir n et e en u et v
		        v_courant = e_courant*cos(deg2rad(cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1)))) - (n_courant*sin(deg2rad(cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1)))));
		        u_courant = -(n_courant*cos(deg2rad(cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1))))) - (e_courant*sin(deg2rad(cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1)))));

		        % Euler forward : interpolation de position
		        cell_data_arbres{3,2}(z+1,ID_arbres_drift_ocean(i,1)) = cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)) + (u_courant*86400);
		        cell_data_arbres{4,2}(z+1,ID_arbres_drift_ocean(i,1)) = cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)) + (v_courant*86400);
		   
		        % Iteration
		        z = z+1;
            end
        end
        [cell_data_arbres{2,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{1,2}(z,ID_arbres_drift_ocean(i,1))] = xytone_V2(cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)));
        cell_data_arbres{5,2}(z,ID_arbres_drift_ocean(i,1)) = interpol_conc_glace_V3(cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)),data_for_interpol, 'linear');
        cell_data_arbres{6,2}(z,ID_arbres_drift_ocean(i,1)) = distance_cote_V2(cote_ORCA, cell_data_arbres{3,2}(z,ID_arbres_drift_ocean(i,1)), cell_data_arbres{4,2}(z,ID_arbres_drift_ocean(i,1)), 1);
        it = i/length(arbres_a_drifter);        
    end
    % Apres l'interpoaltion de la trajectoire via les courants oceaniques, supprimer les jours ou l'arbre est trop pres de la cote (seuil cote);
    cell_data_arbres = filtre_cote(cell_data_arbres,seuil_cote,nb_arbres*2);
    
    disp('Working on formate data and export outputs ...')

    % Une fois que le modele a tourne, diviser le cell array en matrice pour chaque arbre
    % Une arbre est valide sous 2 conditions :
            % 1/ distance entre la position initiale et finale est d'au moins "seuil_exclusion"
            % 2/ loin du point siberien ou les arbres canadiens arrivent
    
    cell_arbres_actifs = cell(1,1);
    % Une ligne = 1 arbre actif
    % La matrice se compose comme suit :
    % 1/ LON
    % 2/ LAT
    % 3/ Concentration glace
    % 4/ Distance cote 
    actif_trees = 0;
    for i = 1:nb_arbres*2;
        % Definir la derniere donnee pour chaque arbre (numero de ligne)
        no_nan = find(isnan(cell_data_arbres{1,2}(:,i)) == 0);
        last_data = max(no_nan);
                
        % Distance entre le point source et le point arrive
        distance_x = cell_data_arbres{3,2}(last_data,i) - cell_data_arbres{3,2}(1,i);
        distance_y = cell_data_arbres{4,2}(last_data,i) - cell_data_arbres{4,2}(1,i);
        distance_abs = sqrt(distance_x^2 + distance_y^2);

        % Une condition supplemntaire !! : Vu que bcp d'arbres se digire au point siberien le plus
        % proche du canada, il faut les retirer 1(ils ne sont pas vrmt actifs..)
        coord_lon_lat_siberie_bug = [180 70];
        [x_siberie_bug, y_siberie_bug] = netoxy(coord_lon_lat_siberie_bug(1), coord_lon_lat_siberie_bug(2));
        distance_x_siberie = cell_data_arbres{3,2}(last_data,i) - x_siberie_bug;
        distance_y_siberie = cell_data_arbres{4,2}(last_data,i) - y_siberie_bug;
        distance_siberie   = sqrt(distance_x_siberie^2 + distance_y_siberie^2);

        if ((distance_abs > distance_exclusion*10^3) && (distance_siberie > 800*10^3));
           actif_trees = actif_trees + 1;
           cell_arbres_actifs{actif_trees,1}(:,1) = cell_data_arbres{1,2}(1:last_data-1,i);
           cell_arbres_actifs{actif_trees,1}(:,2) = cell_data_arbres{2,2}(1:last_data-1,i);
           cell_arbres_actifs{actif_trees,1}(:,3) = cell_data_arbres{5,2}(1:last_data-1,i);
           cell_arbres_actifs{actif_trees,1}(:,4) = cell_data_arbres{6,2}(1:last_data-1,i);
        end
    end 
    
   % Sauvegarder la matrice finale
   if isdir(name_EXP) == 0;
      mkdir(name_EXP);
   end
   if isdir(sprintf('%s/outputs_tmp',name_EXP)) == 0;
      mkdir(sprintf('%s/outputs_tmp',name_EXP));
   end

   save(sprintf('%s/outputs_tmp/TMP_cell_wave_%02d.mat',name_EXP,s),'cell_arbres_actifs')
   % Constuire la matrice qui contient le resume de la simulation
   data_recap = nan(actif_trees,16);   
   % Faire un recap global

                    % 1/num arbre
                    % 2/lon depart
                    % 3/lat depart
                    % 4/lon arrive
                    % 5/lat arrive
                    % 6/duree de la derive
                    % 7/annee des donnees utilisees
                    % 8/duree initiale de la simulation (annees) 
                    % 9/seuil glace
                    % 10/seuil cote
                    % 11/jour de lance des arbres dans l'annee (jour julien)
                    % 12/limite inferieure de la zone de lancement (km)
                    % 13/limite superieure de la zone de lancement (km)
                    % 14/duree de la derive avec les vitesses de la circulation des oceans
                    % 15/distance par rapport a la cote la plus (arrive)
                    % 16/lieu (1=siberie; 2=canada)
   for i = 1:actif_trees;
       % Calculer la duree reelle de la simulation
       [duree_simu_reelle, n_col] = size(cell_arbres_actifs{i,1});
       % Remplir le recap 
       data_recap(i,1)  = i;
       data_recap(i,2)  = cell_arbres_actifs{i,1}(1,1);
       data_recap(i,3)  = cell_arbres_actifs{i,1}(1,2);
       data_recap(i,4)  = cell_arbres_actifs{i,1}(duree_simu_reelle,1);
       data_recap(i,5)  = cell_arbres_actifs{i,1}(duree_simu_reelle,2);
       data_recap(i,6)  = duree_simu_reelle;
       data_recap(i,7)  = 1;
       data_recap(i,8)  = duree_simul;
       data_recap(i,9)  = seuil_glace*100;
       data_recap(i,10) = seuil_cote/100;
       data_recap(i,11) = jour_debut;
       data_recap(i,12) = zone_lancement(1);
       data_recap(i,13) = zone_lancement(2);
       data_recap(i,14) = duree_derive_conc;
       data_recap(i,15) = cell_arbres_actifs{i,1}(duree_simu_reelle,4);
       if ((data_recap(i,2) > 0) & (data_recap(i,2) < 182));
          data_recap(i,16) = 1; % 1 pour SIBERIE
       else
          data_recap(i,16) = 2; % 2 pour CANADA
       end
   end

   % Exporter la matrice
   save(sprintf('%s/outputs_tmp/TMP_recap_wave_%02d.mat',name_EXP,s),'data_recap')
   if isdir(sprintf('%s/FIGURES',name_EXP)) == 0;
       mkdir(sprintf('%s/FIGURES',name_EXP))
   end
   if strcmp(PLOT_TRAJ_O,'on')
   	   for i = 1:actif_trees
		   % Plot des trajectoires
		   figure('Units', 'centimeters','Position', [0 0 25 18],'Visible',SHOW_FIG);
		   hold on
		   m_proj('Azimuthal Equal-Area','lat',90,'long',-90,'radius',30)
		   m_coast('patch',[.5 .5 .5],'edgecolor','none')
		   m_plot(cell_arbres_actifs{i,1}(:,1), cell_arbres_actifs{i,1}(:,2),'r')
		   m_grid_17('linestyle','-','Ytick',[75 86], 'xaxisloc', 'bottom');
		   hold off
		   set(gcf, 'PaperPositionMode', 'auto');
		   if data_recap(i,16) == 1;
		      reg_dep_ID = 'siberia';
		   else
		   	  reg_dep_ID = 'canada';
		   end
		   print('-dpng',sprintf('%s/FIGURES/traj_%s_wave_%d_%04d.png',name_EXP,reg_dep_ID,s,i))
		   close
	   end
   end
    % Afficher le nombre d'arbres actifs
    %fprintf('Actifs trees : %d trees \n',actif_trees)
end

%----------------------------
%
% At the end of waves simulation, export in proper way -> NETCDF
%
%----------------------------

keep seuil_cote seuil_cote duree_derive_conc nb_waves name_EXP
close all
clc

para_coast = seuil_cote/1000;

nb_days = (30*360)+180;
start_days         = sort(randi([90 270],1,nb_waves));

% Merge wave outputs
for s = 1:nb_waves;
    % Load trajectories
    load(sprintf('%s/outputs_tmp/TMP_cell_wave_%02d.mat',name_EXP,s));
    % Load recap
    load(sprintf('%s/outputs_tmp/TMP_recap_wave_%02d.mat',name_EXP,s))

    if isempty(cell_arbres_actifs{1}) == 0

      LON_D_tmp   = nan(length(cell_arbres_actifs),nb_days);
      LAT_D_tmp   = nan(length(cell_arbres_actifs),nb_days);
      ICE_D_tmp   = nan(length(cell_arbres_actifs),nb_days);
      COAST_D_tmp = nan(length(cell_arbres_actifs),nb_days);

      LON_DEPART_tmp   = data_recap(:,2);
      LAT_DEPART_tmp   = data_recap(:,3);
      LON_FINISH_tmp   = data_recap(:,4);
      LAT_FINISH_tmp   = data_recap(:,5);
      ORIGN_D_tmp      = data_recap(:,16);
      COAST_FINISH_tmp = data_recap(:,15);

      for d = 1:length(cell_arbres_actifs);
          LON_D_tmp(d,1:size(cell_arbres_actifs{d},1))   = cell_arbres_actifs{d}(:,1);
          LAT_D_tmp(d,1:size(cell_arbres_actifs{d},1))   = cell_arbres_actifs{d}(:,2); 
          ICE_D_tmp(d,1:size(cell_arbres_actifs{d},1))   = cell_arbres_actifs{d}(:,3);
          COAST_D_tmp(d,1:size(cell_arbres_actifs{d},1)) = cell_arbres_actifs{d}(:,4); 
      end
      if s == 1;
         LON_D   = LON_D_tmp;
         LAT_D   = LAT_D_tmp;
         ICE_D   = ICE_D_tmp;
         COAST_D = COAST_D_tmp;

         LON_DEPART   = LON_DEPART_tmp;
         LAT_DEPART   = LAT_DEPART_tmp;
         LON_FINISH   = LON_FINISH_tmp;
         LAT_FINISH   = LAT_FINISH_tmp;
         ORIGN_D      = ORIGN_D_tmp;
         COAST_FINISH = COAST_FINISH_tmp;

      else
         LON_D   = [LON_D;LON_D_tmp];
         LAT_D   = [LAT_D;LAT_D_tmp];
         ICE_D   = [ICE_D;ICE_D_tmp];
         COAST_D = [COAST_D;COAST_D_tmp];

         LON_DEPART   = [LON_DEPART;LON_DEPART_tmp];
         LAT_DEPART   = [LAT_DEPART;LAT_DEPART_tmp];
         LON_FINISH   = [LON_FINISH;LON_FINISH_tmp];
         LAT_FINISH   = [LAT_FINISH;LAT_FINISH_tmp];
         ORIGN_D      = [ORIGN_D;ORIGN_D_tmp];
         COAST_FINISH = [COAST_FINISH;COAST_FINISH_tmp];
       end
    end
end

if exist('LON_FINISH') ~= 0

    % Analyse des simulations : 
    % 1/ refaire les plots avec diffenciation arbres echouees et non
    % 2/ analyse de la repartition moyenne

    %*************************************%
    %%      REFAIRE LES PLOTS RECAP      %%
    %*************************************%

    % Ici, on ne refait que les plots recap des sorties du modele. Donc on
    % modifie la position des arbres echoues pour qu'elle coincide avec la
    % cote. On exporte egalement les nouvelles tables (modifies)

    % Charger la matrice des couleurs des cotes
    h = colormap(cbrewer('seq','YlGn',100));
    close

    % Charger les coordonnes geo des cotes
    load('data/coord_zone.mat')

    % Charger les traits de cote avec l'ID (3e colonne). Transformer les lon de
    % 0 a 360 et garder moins de points de cotes.
    % 1/ Groenland Sud
    % 2/ Groenland Nord
    % 3/ Islande
    % 4/ Ellesmere
    % 5/ Svalbard
    data_cote_id = load('data/data_cote_ID.mat');
    data_cote_id = data_cote_id.data;
    data_cote_id(data_cote_id(:,1) < 0) = data_cote_id(data_cote_id(:,1) < 0) + 360;
    data_cote_id = data_cote_id(round(linspace(1,80000,1000)),:);

    % Remove driftwood arriving in East Siberia
    lignes_non_sib = find((LON_FINISH > 180) | (LON_FINISH < 35));

    LON_DEPART   = LON_DEPART(lignes_non_sib);
    LAT_DEPART   = LAT_DEPART(lignes_non_sib);
    LON_FINISH   = LON_FINISH(lignes_non_sib);
    LAT_FINISH   = LAT_FINISH(lignes_non_sib);
    LON_D        = LON_D(lignes_non_sib,:);
    LAT_D        = LAT_D(lignes_non_sib,:);
    ICE_D        = ICE_D(lignes_non_sib,:);
    COAST_D      = COAST_D(lignes_non_sib,:);
    COAST_FINISH = COAST_FINISH(lignes_non_sib,:);

    % Remove driftwood arriving in West Canada
    lignes_non_sud_canada = find((LON_FINISH > 270) | (LON_FINISH < 90));
    LON_DEPART   = LON_DEPART(lignes_non_sud_canada);
    LAT_DEPART   = LAT_DEPART(lignes_non_sud_canada);
    LON_FINISH   = LON_FINISH(lignes_non_sud_canada);
    LAT_FINISH   = LAT_FINISH(lignes_non_sud_canada);
    LON_D        = LON_D(lignes_non_sud_canada,:);
    LAT_D        = LAT_D(lignes_non_sud_canada,:);
    ICE_D        = ICE_D(lignes_non_sud_canada,:);
    COAST_D      = COAST_D(lignes_non_sud_canada,:);
    COAST_FINISH = COAST_FINISH(lignes_non_sud_canada,:);
            
    % Compute distance from the nearest coast to assign the coordinate coast
    % to arrival point of driftwood
    supp_dist = 25;
    REG_FINAL_ECH = nan(size(COAST_D,1),1);
    n_arbres = size(LON_FINISH,1);
    for n = 1:n_arbres;
        if COAST_FINISH(n) < ((para_coast + supp_dist)*1000);
            distance_coast = nan(length(data_cote_id),1);
            for j = 1:length(data_cote_id);
                distance_coast(j) = distance_2_pts(LON_FINISH(n),LAT_FINISH(n),data_cote_id(j,1),data_cote_id(j,2));
            end
            % Trouver le minimum
            min_distance = min(distance_coast);
            ligne_min_distance = find(min_distance == distance_coast);

            % Remplacer les coordonnees. Si la cote la plus proche est a
            % plus de 300 km, considere l'arbre comme non echoue
            if (min_distance < 300)
                LON_FINISH(n) = data_cote_id(ligne_min_distance,1);
                LAT_FINISH(n) = data_cote_id(ligne_min_distance,2);
            end           
            % Place of the ending up
            REG_FINAL_ECH(n) = data_cote_id(ligne_min_distance,3);
        else
            REG_FINAL_ECH(n) = 0;
        end
    end
            
    %***************************************************************************************************%
    %%      NOMBRES D'ARBRES POUR CHAQUE REGION, POUR CHAQUE ANNEE ET PROPORTION D'ARBRES ECHOUES      %%
    %***************************************************************************************************% 

    data = nan(size(REG_FINAL_ECH,1),17);
    data(:,2)  = LON_DEPART;
    data(:,3)  = LAT_DEPART;
    data(:,4)  = LON_FINISH;
    data(:,5)  = LAT_FINISH;
    driftwood_canada = find(LON_DEPART>180);
    data(driftwood_canada,16)  = 2;
    data(data(:,16) ~= 2,16) = 1; 
    data(:,17) = REG_FINAL_ECH;

    % Calculer le nombre d'arbres total actif pour chaque region
    data_actifs_sib = data(data(:,16) == 1,:);
    data_actifs_canada = data(data(:,16) == 2,:);
    tot_canada  = length(driftwood_canada);
    tot_siberia = length(REG_FINAL_ECH) - tot_canada;

    % Calculer les propotions
    lignes_arbres_echoues_sib    = find(data_actifs_sib(:,17) ~= 0);
    lignes_arbres_echoues_canada = find(data_actifs_canada(:,17) ~= 0);
    driftwood_siberia_ending_up_prop = 100*length(lignes_arbres_echoues_sib)/lignes_arbres_echoues_sib;
    driftwood_canada_ending_up_prop  = 100*length(lignes_arbres_echoues_canada)/lignes_arbres_echoues_canada;

    %************************************************************%
    %%      NOMBRES D'ARBRES POUR CHAQUE REGION D'ECHOUAGE      %%
    %************************************************************%        

    % Boucler sur les 2 regions (Siberie et Canada)
    for r = 1:2;           
        % Ne garder que les arbres echoues et ceux provenant de la region en question
        data_echoues = data(((data(:,17) ~=0) & (data(:,16) == r)),:);
        % Creer la matrice vide ou l'on stocke la proportion des arbres
        prop_arbres = nan(5,2);            
        % Boucler sur l'ensemble des zones d'arrivages d'arbres
        k = 1;
        for z = [5 2 1 4 3];
            prop_arbres(k,1) = z;
            % Trouver les lignes du fichiers pour la region d'origine et pour la zone d'echouage
            lignes_r_z = find((data_echoues(:,16) == r) & (data_echoues(:,17) == z));
            % Stocker la proportion
            if isempty(lignes_r_z) == 1;
               prop_arbres(k,2) = 0;
            else
               prop_arbres(k,2) = 100*length(lignes_r_z)/length(data_echoues);
               % Round (ne pas avoir des proportions non entieres
               prop_arbres(k,2) = round(prop_arbres(k,2));
            end
            % Compteur
            k = k+1;
        end
        if r == 1;
            prop_arbres_sib = prop_arbres;
        else
            prop_arbres_canada = prop_arbres;
        end
    end

    % Plot du recap de la simulation pour le canada
    figure('Units', 'centimeters','Position', [0 0 10 10],'Visible',SHOW_FIG);
    hold on
    m_proj('Azimuthal Equal-Area','lat',90,'long',-90,'radius',33)
    % Plolter la proportion d'arbres pour chaque cote
    caxis([0 100]);    
    colormap(cbrewer('seq','YlGn',100));
    shading interp
    for c = 1:5;
        if prop_arbres_canada(c,2) == 0;
            % Ne rien plotter
        else
            m_plot(coord_zone{c,2}(:,1), coord_zone{c,2}(:,2),'color',[h(prop_arbres_canada(c,2),1) h(prop_arbres_canada(c,2),2) h(prop_arbres_canada(c,2),3)],'LineStyle','-','LineWidth',8);
        end
    end
    m_coast('patch',[.5 .5 .5],'edgecolor','none')
    % Departs
    m_plot(data(data(:,16) == 2  & data(:,17) ~= 0 ,2) ,data(data(:,16) == 2  & data(:,17) ~= 0,3),'k','LineStyle','none','Marker','.','MarkerSize',8)   
    % Arrives    
    m_plot(data(((data(:,16) == 2) & (data(:,17) ~= 0)),4) ,data(((data(:,16) == 2) & (data(:,17) ~= 0)),5),'r','LineStyle','none','Marker','.','MarkerSize',8)
    %m_plot(data(((data(:,16) == 2) & (data(:,17) == 0)),4) ,data(((data(:,16) == 2) & (data(:,17) == 0)),5),'b','LineStyle','none','Marker','.','MarkerSize',8)
    m_grid_17('linestyle','-','Ytick',[75 75], 'xaxisloc', 'bottom');
    axis equal
    title('CANADA')
    hold off
    set(gcf,'color','w')
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc',sprintf('%s/RECAP_CANADA.eps',name_EXP))

    % Plot du recap de la simulation pour la siberie
    figure('Units', 'centimeters','Position', [0 0 10 10],'Visible',SHOW_FIG);
    hold on
    m_proj('Azimuthal Equal-Area','lat',90,'long',-90,'radius',33)
    % Plolter la proportion d'arbres pour chaque cote
    caxis([0 100]);    
    colormap(cbrewer('seq','YlGn',100));
    shading interp
    for c = 1:5;
        if prop_arbres_sib(c,2) == 0;
            % Ne rien plotter
        else
            m_plot(coord_zone{c,2}(:,1), coord_zone{c,2}(:,2),'color',[h(prop_arbres_sib(c,2),1) h(prop_arbres_sib(c,2),2) h(prop_arbres_sib(c,2),3)],'LineStyle','-','LineWidth',8);
        end
    end
    m_coast('patch',[.5 .5 .5],'edgecolor','none')
    % Departs
    m_plot(data(data(:,16) == 1  & data(:,17) ~= 0,2) ,data(data(:,16) == 1  & data(:,17) ~= 0,3),'k','LineStyle','none','Marker','.','MarkerSize',8)
    %m_plot(data(data(:,16) == 1,2) ,data(data(:,16) == 1,3),'k','LineStyle','none','Marker','.','MarkerSize',8)  
    % Arrives    
    m_plot(data(((data(:,16) == 1) & (data(:,17) ~= 0)),4) ,data(((data(:,16) == 1) & (data(:,17) ~= 0)),5),'r','LineStyle','none','Marker','.','MarkerSize',8)
    %m_plot(data(((data(:,16) == 1) & (data(:,17) == 0)),4) ,data(((data(:,16) == 1) & (data(:,17) == 0)),5),'b','LineStyle','none','Marker','.','MarkerSize',8)
    m_grid_17('linestyle','-','Ytick',[75 75], 'xaxisloc', 'bottom');
    title('SIBERIA')
    axis equal
    hold off
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'color','w')
    print('-depsc',sprintf('%s/RECAP_SIBERIA.eps',name_EXP))       

    path_out = sprintf('%s/%s_OUTPUT.nc',name_EXP,name_EXP);

    nco = netcdf.create(path_out,'CLOBBER');

    dim_id1      = netcdf.defDim(nco,'nb_driftwood',size(LON_D,1));
    dim_id2      = netcdf.defDim(nco,'nb_days',size(LON_D,2));
    dim_id3      = netcdf.defDim(nco,'nb_waves',size(start_days,2));
    dim_region_o = netcdf.defDim(nco,'source_area',2);
    dim_region_f = netcdf.defDim(nco,'arrival_area',5);
    type = 'float';

    var0 = netcdf.defVar(nco,'START_DAYS',type,[dim_id3]);
    netcdf.putAtt(nco,var0,'long_name','Start days of each wave');
    netcdf.putAtt(nco,var0,'units',' ')

    var1 = netcdf.defVar(nco,'LON_DEPART',type,[dim_id1]);
    netcdf.putAtt(nco,var1,'long_name','Depart longitude of each driftwood');
    netcdf.putAtt(nco,var1,'units','degrees_east')
    netcdf.putAtt(nco,var1,'missing_value','Nan');

    var2 = netcdf.defVar(nco,'LAT_DEPART',type,[dim_id1]);
    netcdf.putAtt(nco,var2,'long_name','Depart latitude of each driftwood');
    netcdf.putAtt(nco,var2,'units','degrees_north')
    netcdf.putAtt(nco,var2,'missing_value','Nan');

    var3 = netcdf.defVar(nco,'LON_ARRIVAL',type,[dim_id1]);
    netcdf.putAtt(nco,var3,'long_name','Arrival longitude of each driftwood');
    netcdf.putAtt(nco,var3,'units','degrees_east')
    netcdf.putAtt(nco,var3,'missing_value','Nan');

    var4 = netcdf.defVar(nco,'LAT_ARRIVAL',type,[dim_id1]);
    netcdf.putAtt(nco,var4,'long_name','Arrival latitude of each driftwood');
    netcdf.putAtt(nco,var4,'units','degrees_north')
    netcdf.putAtt(nco,var4,'missing_value','Nan');

    var5 = netcdf.defVar(nco,'DIST_COAST_ARRIVAL',type,[dim_id1]);
    netcdf.putAtt(nco,var5,'long_name','Distance from the nearest coast at arrival point');
    netcdf.putAtt(nco,var5,'units','m')
    netcdf.putAtt(nco,var5,'missing_value','Nan');

    var10 = netcdf.defVar(nco,'LON_TIME_SERIE',type,[dim_id2 dim_id1]);
    netcdf.putAtt(nco,var10,'long_name','Longitude of each driftwood at each time step');
    netcdf.putAtt(nco,var10,'units','degrees_east')
    netcdf.putAtt(nco,var10,'missing_value','Nan');

    var20 = netcdf.defVar(nco,'LAT_TIME_SERIE',type,[dim_id2 dim_id1]);
    netcdf.putAtt(nco,var20,'long_name','Latitude of each driftwood at each time step');
    netcdf.putAtt(nco,var20,'units','degrees_north')
    netcdf.putAtt(nco,var20,'missing_value','Nan');

    var30 = netcdf.defVar(nco,'ICE_CONC_TIME_SERIE',type,[dim_id2 dim_id1]);
    netcdf.putAtt(nco,var30,'long_name','Sea ice concentration of each driftwood position at each time step');
    netcdf.putAtt(nco,var30,'units','/')
    netcdf.putAtt(nco,var30,'min','0')
    netcdf.putAtt(nco,var30,'max','1')
    netcdf.putAtt(nco,var30,'missing_value','Nan');

    var40 = netcdf.defVar(nco,'COAST_DIST_TIME_SERIE',type,[dim_id2 dim_id1]);
    netcdf.putAtt(nco,var40,'long_name','Minimal distance from the neartest coast of each driftwood position at each time step');
    netcdf.putAtt(nco,var40,'units','m')
    netcdf.putAtt(nco,var40,'missing_value','Nan');

    var50 = netcdf.defVar(nco,'prop_canada_driftwood',type,[dim_region_o dim_region_f]);
    netcdf.putAtt(nco,var50,'long_name','Proportion of Canadian driftwood for each arrival region. 1 = South Greenland; 2 = North Greenland; 3 = Iceland; 4 = Svalbard; 5 = Ellesmere Island');
    netcdf.putAtt(nco,var50,'units','%')
    netcdf.putAtt(nco,var50,'missing_value','Nan');

    var60 = netcdf.defVar(nco,'prop_siberia_driftwood',type,[dim_region_o dim_region_f]);
    netcdf.putAtt(nco,var60,'long_name','Proportion of Siberian driftwood for each arrival region. 1 = South Greenland; 2 = North Greenland; 3 = Iceland; 4 = Svalbard; 5 = Ellesmere Island');
    netcdf.putAtt(nco,var60,'units','%')
    netcdf.putAtt(nco,var60,'missing_value','Nan');

    % Global attributes
    varid_G = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(nco,varid_G,'decription','Driftwood Transport Model (DTM)')
    netcdf.putAtt(nco,varid_G,'Parameters',sprintf('coast_thr=%d km; ice_thr=%0.2f; ice_time=%d days',seuil_cote,seuil_glace,duree_derive_conc))
    netcdf.putAtt(nco,varid_G,'contact','Dalaiden Quentin (quentin.dalaiden@uclouvain.be)')
    netcdf.putAtt(nco,varid_G,'insitute','Universite catholique de Louvain (UCL)')
    netcdf.putAtt(nco,varid_G,'create_file',date)

    netcdf.endDef(nco)

    start_days = permute(start_days,[2 1]);
    netcdf.putVar(nco,var0,start_days)
    LON_DEPART = permute(LON_DEPART,[2 1]);
    netcdf.putVar(nco,var1,LON_DEPART)
    LAT_DEPART = permute(LAT_DEPART,[2 1]);
    netcdf.putVar(nco,var2,LAT_DEPART)
    LON_FINISH = permute(LON_FINISH,[2 1]);
    netcdf.putVar(nco,var3,LON_FINISH)
    LAT_FINISH = permute(LAT_FINISH,[2 1]);
    netcdf.putVar(nco,var4,LAT_FINISH)
    COAST_FINISH = permute(COAST_FINISH,[2 1]);
    netcdf.putVar(nco,var5,COAST_FINISH)


    LON_D = permute(LON_D,[2 1]);
    netcdf.putVar(nco,var10,LON_D)
    LAT_D = permute(LAT_D,[2 1]);
    netcdf.putVar(nco,var20,LAT_D)
    ICE_D = permute(ICE_D,[2 1]);
    netcdf.putVar(nco,var30,ICE_D)
    COAST_D = permute(COAST_D,[2 1]);
    netcdf.putVar(nco,var30,COAST_D)
    prop_arbres_canada = permute(sortrows(prop_arbres_canada,1),[2 1]);
    netcdf.putVar(nco,var50,prop_arbres_canada)
    prop_arbres_sib = permute(sortrows(prop_arbres_sib,1),[2 1]);
    netcdf.putVar(nco,var60,prop_arbres_sib)

    netcdf.close(nco)

else
    fid = fopen(sprintf('%s/%s_OUTPUT.txt',name_EXP,name_EXP),'w');
    fprintf(fid,'no driftwood for this experience');
    fclose(fid)
end

end
