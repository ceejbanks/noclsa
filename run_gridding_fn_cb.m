% function that does the meat in the sandwich of grid_data_by_region

function        run_gridding_fn_cb(fn, lons_grid,lats_grid,  ...
    SR_distance,SR_time,E,ef_distance,ef_time,days_use,vars_to_use,wrap_lon)
[fn_path,fn_out,~] = fileparts(fn);
if ~strcmp(fn_path(end),'/')
    fn_path =[fn_path '/'];
end
fn_out = [fn_path 'GRIDDED_' fn_out '.mat'];
% load datafile
load(fn)
% now need to reduce size of grids so only have in-range data if not done
% before
n_files = length(days_use);
% initialise
median_ssha = NaN([size(lats_grid) n_files]) ;
mean_ssha = median_ssha ;
weighted_mean_ssha = median_ssha ;
nobs = median_ssha;
sd_out = median_ssha ;
weighted_sd = median_ssha ;
weighted_median = median_ssha ;
[YY,MM,DD] =  datevec(tutc) ;
date_group = datenum(YY,MM,DD) ;
if wrap_lon
    lon = wrapTo360(lon) ;
end
parpool(24)
parfor ntime = (SR_time + 1) :( length(days_use) - SR_time)
    these_time =   (date_group > ( floor(days_use(ntime)) - SR_time) & date_group <= ( ceil(days_use(ntime)) + SR_time)) ;
    if  sum(these_time) >0
        [weighted_mean_ssha(:,:,ntime) , median_ssha(:,:,ntime) , ...
            mean_ssha(:,:,ntime),   nobs(:,:,ntime), sd_out(:,:,ntime),...
            weighted_sd(:,:,ntime),weighted_median(:,:,ntime)] = ...
            loopy_calcs_RROD_BASIN(E, SR_time ,  ef_time ,...
            date_group(these_time) , lats_grid , lons_grid , ...
            lat(these_time),lon(these_time) , ssha(these_time),...
            days_use(ntime) ,SR_distance,ef_distance ) ;
    end
end; clear ntime

for nvars = 1:length(vars_to_use)
    eval(['clear ' vars_to_use{nvars}])
end; clear nvars
poolobj = gcp('nocreate');
delete(poolobj);
clear poolobj

clear YY MM DD date_group E ef_distance ef_time fn fn_path n n_files ...
    date_out  pass_id these_time SR_distance SR_time use_ocean_basin ans vars_to_use
save(fn_out,'-v7.3')











