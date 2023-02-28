function     [weighted_mean_ssha , median_ssha , mean_ssha,...
    nobs , sd_out , weighted_sd , weighted_median] =  ...
    loopy_calcs_RROD_BASIN(E, ...
    SR_time ,  ef_time , date_group, lats_grid2,...
    lons_grid ,lat , lon , ssha_temp , this_date ,SR_distance,ef_distance)
% initialise
median_ssha = NaN(size(lons_grid) ) ;
mean_ssha = median_ssha ;
weighted_mean_ssha = median_ssha ;
nobs = median_ssha;
sd_out = median_ssha ;
weighted_sd  = median_ssha ;
[nrow_max , ncol_max  ] = size(lons_grid) ;

datelen_temp = abs(date_group - this_date) ;

parfor nrow = 1:nrow_max ; % 621
    % ASSUMES is regaulr grid which it normally is
    lats_grid = lats_grid2(nrow,1) ;
    for ncol = 1:ncol_max ; % 773
        if ~isnan(SR_distance(nrow , ncol)) ;
            these_cols  = ...
                (lat <= ceil(lats_grid + (SR_distance(nrow,ncol) ./ 111)) ) & ...
                (lat >= floor(lats_grid - (SR_distance(nrow,ncol) ./ 111)) ) & ...
                (lon <= ceil(lons_grid(nrow,ncol) + (SR_distance(nrow,ncol) ./ 111)) ) & ...
                (lon >= floor(lons_grid(nrow,ncol) - (SR_distance(nrow,ncol) ./ 111)) ) & ...
                ((datelen_temp / SR_time) .^ 2) < 1;
            
            if sum(these_cols) > 0
                datelen = datelen_temp(these_cols) ;
                arclen = distance(lats_grid,lons_grid(nrow,ncol),lat(these_cols),lon(these_cols) ,E )  ./ 1000 ;
                %                 arclen = arclen   ;
                % find within SR_distance
                these_dist = ((arclen / SR_distance(nrow,ncol)) .^ 2) < 1 ;
                % split by number of obs as do not need weighting for where
                % just one measurement
                
                if sum(these_dist) > 0
                    these_dist = ((arclen / SR_distance(nrow,ncol)) .^ 2) + ...
                        ((datelen / SR_time) .^ 2) < 1 ;
                    if sum(these_dist) > 0
                        ssha  = ssha_temp(these_cols) ;
                        ssha =  ssha(these_dist) ;
                        nobs(nrow,ncol) = sum(~isnan(ssha)) ;
                        if nobs(nrow,ncol) == 1
                            try
                                median_ssha(nrow,ncol) = ssha ;
                                mean_ssha(nrow,ncol) = ssha ;
                                sd_out(nrow,ncol) = 0 ;
                                weighted_mean_ssha(nrow,ncol) = 0;
                                weighted_sd(nrow,ncol) = 0;
                                weighted_median(nrow,ncol) = 0;
                                ssha = [];
                                
                            catch
                                median_ssha(nrow,ncol) = NaN ;
                            end
                            
                        elseif nobs(nrow,ncol) > 1
                            % calculate weights
                            W = return_w_global(arclen(these_dist) ,datelen(these_dist) ,ef_distance(nrow,ncol) , ef_time) ;
                            median_ssha(nrow,ncol) = nanmedian(ssha) ;
                            mean_ssha(nrow,ncol) = nanmean(ssha) ;
                            sd_out(nrow,ncol) = nanstd(ssha) ;
                            [weighted_mean_ssha(nrow,ncol), weighted_sd(nrow,ncol) ,...
                                weighted_median(nrow,ncol)] = ...
                                weighted_mean_etc_fn(ssha,W) ;
                            w_time =[]; w_distance =[];arclen=[]; these_cols = [];%#ok<NASGU>
                            ssha = [];
                        end ; % end if on NaN lon
                    end;
                    these_dist = [] ;arclen = [] ; datelen = []; %clear sum_these_dist
                end;
            end;
        end
    end ; %clear ncol
end; %clear nrow
