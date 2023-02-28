%  script to look at building DAILY L4 gridded product for Cryosat-2
clear;close all; clc
save_to = data_location';
load for_gridding.mat
FWHM_time = 15 ;
ros_rad_fact = 2 ; % by default is 2, if > 5 is fixed value
% grid size in MINUTES
grid_size_min = 15 ;
E = referenceEllipsoid('wgs84') ;
SR_time = ceil((FWHM_time /2 ) * 3) ;
ef_time=(FWHM_time./2)./sqrt(log(2));
days_use = start_date:end_date;
dates_use = days_use; 
num_days = length(dates_use);
save dates_used_file2020.mat days_use num_days
clear num_days dates_use
% runs over (over lapping) regions
for nregion =  1:8
    disp(['running ' num2str(nregion)])
    datestr(now)
    fn = [file_inut_folder 'global' num2str(nregion)];
        if ismember(nregion,[1 4 5 8])
        wrap_lon = true;
        if ismember(nregion,[1  5 ])
            max_lon_use(nregion) = 270;   min_lon_use(nregion) = 180;
        end
        if ismember(nregion,[4 8 ])
            max_lon_use(nregion) = 180;     min_lon_use(nregion) = 90;
        end
    else
        max_lon_use(nregion) = max_lon(nregion)   ;
        min_lon_use(nregion) = min_lon(nregion);
        wrap_lon = false;
    end
    if ismember(nregion,[1 2 3 4])
        max_lat_use(nregion) = 60 ;        min_lat_use(nregion) = 0 ;
    else
        max_lat_use(nregion) = 0 ;        min_lat_use(nregion) = -60 ;
    end
    % calculate grids
    lons_grid = min_lon_use(nregion):(grid_size_min / 60):(max_lon_use(nregion) -(grid_size_min / 60 )) ;
    nlon = length(lons_grid) ;
    lats_grid = min_lat_use(nregion):(grid_size_min / 60):(max_lat_use(nregion) -(grid_size_min / 60 )) ;
    nlat = length(lats_grid) ;    lons_grid = repmat(lons_grid,[nlat 1]) ;
    lats_grid = repmat(lats_grid',[1 nlon]);
    if ros_rad_fact > 5
        % fixed vlaue for length scale
        ros_rad = lons_grid + NaN ;        ros_rad(:) = ros_rad_fact;
        FWHM_distance =  ros_rad;
    else
        % calculate Rossby radius of deformation
        if ~wrap_lon
            % function to return RRoD from Chelton
            ros_rad =  ross_radius_fn_global(lons_grid , lats_grid) ;
            % see which lons/lats are land (so can be ignored), based on
            % SMOS distance to land product
            use_ocean_basin =  get_dist_from_land_cryo(lons_grid , lats_grid) == 1;
        else
            if ismember(nregion,[1 5 ])
                lon_temp = wrapTo180(lons_grid) ;
                lon_temp(:,1) = -180;
                ros_rad =  ross_radius_fn_global(lon_temp , lats_grid) ;
                % see which lons/lats are land (so can be ignored)
                use_ocean_basin =  get_dist_from_land_cryo(lon_temp , lats_grid) == 1;
                clear lon_temp
            elseif ismember(nregion,[4 8 ])
                % see which lons/lats are land (so can be ignored)
                use_ocean_basin =  get_dist_from_land_cryo(lons_grid , lats_grid) == 1;
                ros_rad =  ross_radius_fn_global(wrapTo180(lons_grid) , lats_grid) ;
            else
                error('grrr')
            end;
        end
        FWHM_distance = ros_rad_fact .* ros_rad ;
    end
    SR_distance = ceil((FWHM_distance ./2 ) .* 3) ;
    ef_distance = (FWHM_distance./2)./sqrt(log(2));
    max_SR_dist = max(SR_distance(:)) ;
    SR_distance(~use_ocean_basin) = NaN;
    run_gridding_fn_cb(fn, lons_grid,lats_grid,  ...
        SR_distance,SR_time,E,ef_distance,ef_time,days_use,vars_to_use,wrap_lon);
    disp(['Have run ' num2str(nregion)])
    datestr(now)
    poolobj = gcp('nocreate');
    delete(poolobj);
    clear poolobj
end


