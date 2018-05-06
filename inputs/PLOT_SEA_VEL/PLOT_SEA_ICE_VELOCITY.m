% Plot Sea ice velocity

clear all
close all
clc

addpath(genpath('~/MATLAB'))

MODELS         = {'HadGEM2-CC','FGOALS-g2'};

for m = 1%:length(MODELS);

	% Load variables
	[LON,LAT,velx] = load_netcdf(sprintf('/ofast/dalaiden/INPUT_DTM_MODEL/MH/%s/usi_%s_CLIMATOLOGY_MH.nc',MODELS{m},MODELS{m}),'LON','LAT','usi');
	[LON,LAT,vely] = load_netcdf(sprintf('/ofast/dalaiden/INPUT_DTM_MODEL/MH/%s/vsi_%s_CLIMATOLOGY_MH.nc',MODELS{m},MODELS{m}),'LON','LAT','vsi');

	% Compute drift
	drift = nan(size(velx));
	for t = 1:size(velx,1);
		drift(t,:,:) = sqrt((velx(t,:,:).^2) + (vely(t,:,:).^2)).*86.4;
	end

	% Interpolation of the vectors
	LON_vecteurs = load(sprintf('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/LON_eq_%d.mat',23)); 
	LON_vecteurs = LON_vecteurs.LON_eq;
	LAT_vecteurs = load(sprintf('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/LAT_eq_%d.mat',23)); 
	LAT_vecteurs = LAT_vecteurs.LAT_eq;
	load('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/masque_vect_2.mat')

	mean_e = squeeze(nanmean(velx));
	mean_n = squeeze(nanmean(vely));

	mean_e_I = interpol_grid(LON,LAT,mean_e,LON_vecteurs,LAT_vecteurs);
	mean_n_I = interpol_grid(LON,LAT,mean_n,LON_vecteurs,LAT_vecteurs);	

	[row_sup_70, col_sup_70] = find((LAT_vecteurs < 70) & (LON_vecteurs > 300));
	mean_e_I(row_sup_70,col_sup_70) = nan;
	mean_n_I(row_sup_70,col_sup_70) = nan;

	% interpolation of the drift
	load('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/LON_eq_51.mat');
	load('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/LAT_eq_51.mat');
	load('/efast/dalaiden/MASTER_THESIS/biais_modele_transport/masque_vecteurs_drift/masque_drift.mat');
	%[r_beug_ellesmere,c_beug_ellesmere] = find((LAT > 80) & (LAT < 85) & (LON(:,2) > 285) & (LON(:,2) < 315));
	mean_drift = squeeze(nanmean(drift));
	mean_drift_I = interpol_grid(LON,LAT,mean_drift,LON_eq,LAT_eq);
	%mean_drift(r_beug_ellesmere,c_beug_ellesmere) = nan;

	% Figure
	figure%('Units', 'centimeters','Position', [0 0 8.5 14.2]);
	hold on
	m_proj('Azimuthal Equal-Area','lat',90,'long',0,'radius',25)
	%m_plot(LON(b == 1),LAT(b == 1),'linestyle','none','marker','.','markersize',8)
	m_pcolor(LON_eq,LAT_eq,mean_drift_I.*masque_drift)
	colormap_MAT = cbrewer('seq','Reds',60);
	colormap_MAT(1,:) = 1;
	colormap(colormap_MAT)
	shading interp
	m_coast('patch',[.5 .5 .5],'edgecolor',[.5 .5 .5])
	m_vec(0.2,LON_vecteurs, LAT_vecteurs, mean_e_I.*masque_vect_2,mean_n_I.*masque_vect_2,'k', 'shaftwidth', .45, 'headangle', 35, 'headlength', 4)
	caxis([0 15])
	colorbar
	m_grid_17('linestyle','--','Ytick',[75 75], 'xaxisloc', 'bottom');
end


