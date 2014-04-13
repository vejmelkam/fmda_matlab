%
% Load station data from a set of CSV files.  Loads the fm10 observations,
% relative humidity, air temperature and accumulated precipitation.
%
% Returns everything in one data structure sd with fields:
%    tdays, fm10, fm10_var, relh, t2, accp
%
% synopsis: sd = load_station_data(station,Ts)
%
%         station - station code
%         Ts - an array of timestamps from WRF
%         R - either not present, [StartDay, EndDay]
%

function sd = load_station_data_for_wrf(station,Ts)

    year = datevec(Ts(1));
    year = num2str(year(1));
    
    % load all csv files
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
        
    % chop to 1..N if desired
    Nt = length(Ts);
    fm10 = nan * ones(Nt,1);
    fm10_var = nan * ones(Nt,1);
    relh = nan * ones(Nt,1);
    t2 = nan * ones(Nt,1);
    accp = nan * ones(Nt,1);
    for i=1:Nt
        t = Ts(i);
        
        [td,si] = min(abs(tdays1 - t));
        if(td < 1/47)
            fm10(i) = fm10_in(si,7);
            fm10_var(i) = fm10_in(si,8);
        end
        [td,si] = min(abs(tdays2 - t));
        if(td < 1/47)
            relh(i) = relh_in(si,7);
        end
        [td,si] = min(abs(tdays3 - t));
        if(td < 1/47)
            t2(i) = t2_in(si,7);
        end
        [td,si] = min(abs(tdays4 - t));
        if(td < 1/47)
            accp(i) = accp_in(si,7);
        end
    end
    
    % remove any zero observations
    fm10_var(fm10==0) = nan;
    fm10(fm10==0) = nan;
    
    % finds maximum and if it is in the dataset more than once, marks it as
    % a no observation (we have reason to believe these are stronly censored)
    % alternative strategy: increase the variance of the observation
    % markedly
    mx = max(fm10(isfinite(fm10)));
    if(~isempty(mx) && (sum(fm10==mx) > 1))
        fm10_var(fm10==mx) = nan;
        fm10(fm10==mx) = nan;
    end

    % convert accumulated precipitation into rain
    rain = [0;max(diff(accp),0) ./ (diff(Ts) * 24)];
    
    sd.tdays = Ts;
    sd.t2 = t2;
    sd.relh = relh;
    sd.fm10 = fm10;
    sd.fm10_var = fm10_var;
    sd.rain = rain;
    sd.accp = accp;
    sd.elevation = info(3);
    sd.lon = info(2);
    sd.lat = info(1);
    
    [ed,ew] = equilibrium_moisture2(relh,t2);
    sd.eds = ed;
    sd.ews = ew;
    
    
    
    
