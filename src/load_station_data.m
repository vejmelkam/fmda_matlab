%
% Load station data from a set of CSV files.  Loads the fm10 observations,
% relative humidity, air temperature and accumulated precipitation.
%
% Returns everything in one data structure sd with fields:
%    tdays, fm10, fm10_var, relh, t2, accp
%
% synopsis: sd = load_station_data(station,year,R)
%
%         station - station code
%         yr - 4digit year string
%         R - either not present, or EndSample or [BeginSample,EndSample]
%

function sd = load_station_data(station,year,R)

    % load all 4 csv files.
    fm10_in = load(['../data/',station,'_fm10_',year,'.csv']);
    relh_in = load(['../data/',station,'_rel_humidity_',year,'.csv']);
    t2_in = load(['../data/',station,'_air_temp_',year,'.csv']);
    accp_in = load(['../data/',station,'_accum_precip_',year,'.csv']);
    info = load(['../data/',station,'_data.csv']);

    % sort by date/time
    fm10_in = sortrows(fm10_in);
    relh_in = sortrows(relh_in);
    t2_in = sortrows(t2_in);
    accp_in = sortrows(accp_in);
    
    % compute equilibria from station data
    tdays1 = datenum(fm10_in(:,1:6));
    tdays2 = datenum(relh_in(:,1:6));
    tdays3 = datenum(t2_in(:,1:6));
    tdays4 = datenum(accp_in(:,1:6));
    
    % find the intersection of all the weather data
%    [tdays12,map1,map2] = intersect(tdays1,tdays2);
    [tdays34,map3,map4] = intersect(tdays3,tdays4);
    [tdays,map2,map34] = intersect(tdays2,tdays34);
    
    % assigns nan if no observation present or the fm10 observation so that
    % fm10/fm10_var is matched with relh,t2,accp
    fm10 = assign_via_tdays(tdays1,fm10_in(:,7),tdays);
    fm10_var = assign_via_tdays(tdays1,fm10_in(:,8),tdays);
    
    % remove any zero observations
    fm10_var(fm10==0) = nan;
    fm10(fm10==0) = nan;
    
    % finds maximum and if it is in the dataset more than once, marks it as
    % a no observation (since these may be strongly censored)
    % alternative strategy: increase the variance of the observation
    % markedly
    mx = max(fm10(isfinite(fm10)));
    if(sum(fm10==mx) > 1)
        fm10_var(fm10==mx) = nan;
        fm10(fm10==mx) = nan;
    end
    
    % process the observations
    relh = relh_in(map2,7);
    t2   = t2_in(map3(map34),7);
    accp = accp_in(map4(map34),7);
    
    % chop to 1..N if desired
    if(nargin == 3)
        if(numel(R) == 0)
            From = 1;
            To = length(tdays);
        elseif(numel(R) == 1)
            From = 1;
            To = R;
        else
            From = R(1);
            To = R(2);
        end
        tdays = tdays(From:To);
        fm10 = fm10(From:To);
        fm10_var = fm10_var(From:To);
        relh = relh(From:To);
        t2 = t2(From:To);
        accp = accp(From:To);
    end
    
    % convert accumulated precipitation into rain
    rain = [0;max(diff(accp),0) ./ (diff(tdays) * 24)];
    
    sd.tdays = tdays;
    sd.t2 = t2;
    sd.relh = relh;
    sd.fm10 = fm10;
    sd.fm10_var = fm10_var;
    sd.rain = rain;
    sd.accp = accp;
    sd.elevation = info(3);
    sd.lon = info(2);
    sd.lat = info(1);
    
    [Ed,Ew] = equilibrium_moisture2(sd.relh,sd.t2);
    sd.ed = Ed;
    sd.ew = Ew;
    sd.stid = station;
    
    
    
