% Cat MH nc

clear all
close all
clc

addpath(genpath('~/MATLAB'))

setenv('LC_ALL','C')

MODELS         = {'HadGEM2-CC','FGOALS-g2'};
VAR_N          = {'vsi','usi','vo','uo','sic'};
VAR_LONG_NAME  = {'Sea ice velocity Y component',...
				  'Sea ice velocity X component',...
				  'Ocean ice velocity Y component',...
				  'Ocean ice velocity X component',...
				  'Sea Ice concentration'};

VAR_UNITS       = {'m/s','m/s','m/s','m/s','0 to 1'};
directory_file = '~/lfast/DOWNLOAD_CMIP5/MH';

for m = 1%1:length(MODELS);
	for v =[3 4]% 1:length(VAR_N);
		fprintf('MODEL : %s \n',MODELS{m});
		fprintf('VAR   : %s \n',VAR_N{v});

		files = rdir(sprintf('%s/%s_*_%s_*.nc',directory_file,VAR_N{v},MODELS{m}));
		for f = 1:length(files);
			fprintf('		%s \n',files(f).name)
			[LAT, LON, VAR] = load_netcdf(files(f).name,'lat','lon',VAR_N{v});
			VAR(VAR >10^10) = nan;
			if strcmp(VAR_N{v},'uo') | strcmp(VAR_N{v},'vo');
				VAR = squeeze(VAR(:,1,:,:));
			elseif strcmp(VAR_N{v},'sic');
				VAR = VAR./100;
			end
			VAR = permute(VAR,[2 3 1]);

			% CAT
			if f == 1;
				CAT_VAR = VAR;
			else
				CAT_VAR = cat(3,CAT_VAR,VAR);
			end
		end

		clear VAR

		%------------------------------
		% Compute the daily climatolgy
		%------------------------------

		% For HadGEM2-CC, remove the first month
		if strcmp(MODELS{m},'HadGEM2-CC')
			if strcmp(VAR_N{v},'uo') | strcmp(VAR_N{v},'vo');
				CAT_VAR = CAT_VAR(:,:,2:end);
			else
				CAT_VAR = CAT_VAR(:,:,31:end);
			end
		end

		% Verify the size
		if strcmp(VAR_N{v},'uo') | strcmp(VAR_N{v},'vo');
			size_NC = size(CAT_VAR,3)/12;
		else
			if strcmp(MODELS{m},'FGOALS-g2');
				size_NC = size(CAT_VAR,3)/365;
			else
				size_NC = size(CAT_VAR,3)/36;
			end
		end

		if mod(size_NC,1) ~= 0;
			disp('problem')
			break
		end

		% Permute the VAR
		CAT_VAR = permute(CAT_VAR,[3 1 2]);

		% compute the climatology
		if strcmp(MODELS{m},'FGOALS-g2') & (strcmp(VAR_N{v},'uo') == 0 | strcmp(VAR_N{v},'vo') == 0);
			CLIM_DAILY = daily_climatology(CAT_VAR,365);
		elseif strcmp(MODELS{m},'HadGEM2-CC') & (strcmp(VAR_N{v},'uo') == 0 & strcmp(VAR_N{v},'vo') == 0);
			CLIM_DAILY = daily_climatology(CAT_VAR,360);
		else
			CLIM_DAILY = monthmean(CAT_VAR);
		end

		% For FGOALS-g2 -> switch to 360 days and not 360 days
		if strcmp(MODELS{m},'FGOALS-g2') & (strcmp(VAR_N{v},'uo') == 0 | strcmp(VAR_N{v},'vo') == 0);
			NEW_CLIM_DAILY = nan(360,size(CLIM_DAILY,2),size(CLIM_DAILY,3));
			for x = 1:size(CLIM_DAILY,2);
				for y = 1:size(CLIM_DAILY,3);
					i_365 = 1:365;
					i_360 = linspace(1,365,360);
					NEW_CLIM_DAILY(:,x,y) = interpn(i_365,squeeze(CLIM_DAILY(:,x,y)),i_360,'linear'); 
				end
			end
			clear CLIM_DAILY
			CLIM_DAILY = NEW_CLIM_DAILY;
			clear NEW_CLIM_DAILY
		end

		% Interpol on NEMO grid
		[LAT_NEMO,LON_NEMO,MASK_NEMO] = load_netcdf('~/lfast/DTM_RAW_INPUTS/mesh_mask_NEMO1.nc','nav_lat','nav_lon','tmask');
		MASK_NEMO = squeeze(MASK_NEMO(1,:,:));
		MASK_NEMO(MASK_NEMO == 0) = nan;

		[LON,LAT] = meshgrid(LON,LAT);

		VAR_INTERP_2 = nan(size(CLIM_DAILY,1),size(LAT_NEMO,1),size(LAT_NEMO,2));
		p = parpool(12)
		parfor m_t = 1:size(CLIM_DAILY,1);
			VAR_INTERP_2(m_t,:,:) = interpol_grid(LON,LAT,squeeze(CLIM_DAILY(m_t,:,:)),LON_NEMO,LAT_NEMO);
			VAR_INTERP_2(m_t,:,:) = inpaintn(squeeze(VAR_INTERP_2(m_t,:,:))).*MASK_NEMO;
		end
		delete(gcp('nocreate'))

		%----------------------
		% Export in netcdf
		%----------------------

		path_out=sprintf('/lfast/dalaiden/INPUT_DTM_MODEL/%s/%s/%s_%s_midHolocene.nc','midHolocene',MODELS{m},VAR_N{v},MODELS{m});
        nco = netcdf.create(path_out,'CLOBBER');
		dim_id1 = netcdf.defDim(nco,'x',size(LAT_NEMO,1));
		dim_id2 = netcdf.defDim(nco,'y',size(LAT_NEMO,2));
		dim_id3 = netcdf.defDim(nco,'time',size(VAR_INTERP_2,1));
		type='float';

		var1 = netcdf.defVar(nco,'LON',type,[dim_id2 dim_id1]);
		netcdf.putAtt(nco,var1,'long_name','longitude');
		netcdf.putAtt(nco,var1,'units','degrees_east');

		var2 = netcdf.defVar(nco,'LAT',type,[dim_id2 dim_id1]);
		netcdf.putAtt(nco,var2,'long_name','latitude');
		netcdf.putAtt(nco,var2,'units','degrees_north');

		var3 = netcdf.defVar(nco,VAR_N{v},type,[dim_id2 dim_id1 dim_id3]);
		netcdf.putAtt(nco,var3,'standart_name',VAR_LONG_NAME{v})
		netcdf.putAtt(nco,var3,'units',VAR_UNITS{v})

		% Global attributes
		varid_G = netcdf.getConstant('GLOBAL');
		netcdf.putAtt(nco,varid_G,'create_file',date)
		netcdf.putAtt(nco,varid_G,'Model',MODELS{m})

		netcdf.putAtt(nco,varid_G,'Author','Quentin Dalaiden')
		netcdf.putAtt(nco,varid_G,'Program_used','MATLAB')

		netcdf.endDef(nco)

		LON_NEMO = permute(LON_NEMO,[2 1]);
		netcdf.putVar(nco,var1,LON_NEMO)
		LAT_NEMO = permute(LAT_NEMO,[2 1]);
		netcdf.putVar(nco,var2,LAT_NEMO)
		netcdf.putVar(nco,var3,permute(VAR_INTERP_2,[3 2 1]))

		netcdf.close(nco)

	end
end

