

%
% This function conducts an experiment evaluating the ability of the TSM to
% transport observation information across space using the Fay-Herriott
% model.  The testing is done via a leave-one-out testing strategy.
%
%  synopsis: out = experiment_tsm_wrf(wrf_file)
%

function experiment_model_wrf(wrf_file)

    % read out data from wrf file
    wrf_t = ncread(wrf_file, 'Times')';
    Nt = size(wrf_t,1);
    psfc = ncread(wrf_file,'PSFC');
    q2 = ncread(wrf_file, 'Q2');
    t2 = ncread(wrf_file,'T2');
    glon = ncread(wrf_file,'XLONG');
    glon = glon(:,:,1);
    glat = ncread(wrf_file,'XLAT');
    glat = glat(:,:,1);
    rain = cat(3,zeros(size(glon)),diff(ncread(wrf_file,'RAINC'),1,3) + diff(ncread(wrf_file,'RAINNC'),1,3));
    [Nlon,Nlat] = size(glon);
    [ed,ew,relh] = equilibrium_moisture(psfc,q2,t2);
    relh = relh / 100;
    
    assert(all(ed(:) > 0));
    assert(all(ew(:) > 0));
    
    if(any(relh(:) < 0))
        warning('fixing WRF relative humidity, found %d negative values', sum(relh(:) < 0));
        relh(relh < 0) = 0;
    end
    if(any(relh(:) > 1))
        warning('fixing WRF relative humidity, found %d negative values', sum(relh(:) > 1));
        relh(relh > 1) = 1;
    end
    
    % parse the dates
    Ts = zeros(Nt,1);
    for i=1:Nt
        ti = wrf_t(i,:);
        Ts(i) = datenum(str2double(ti(1:4)),str2double(ti(6:7)),str2double(ti(9:10)),...
                        str2double(ti(12:13)),str2double(ti(15:16)),str2double(ti(18:19)));
        % convert rain to mm/h units
        if(i>1)
            rain(:,:,i) = rain(:,:,i) / ((Ts(i) - Ts(i-1)) * 24);
        end
    end
    
    % remove any rain that is too light to affect the moisture model
    rain(rain < 0.01) = 0;
    
    % load data from correct year
    year = str2double(ti(1:4));
    years = num2str(year);
    S = dir(['../data/*fm10*',years,'*']);
    N = length(S);

    sds = cell(N,1);
    base_tday = datenum(year,1,1,0,0,0);
    
    % output dir prefix
    ti = wrf_t(1,:);
    out_dir = sprintf('wrf_%d%d%d_model', str2double(ti(1:4)), str2double(ti(6:7)), str2double(ti(9:10)));
    mkdir(out_dir);
    
     % constant model parameters (from fit to first 10 WRF files)
    Tk = [1,10,100]';
    M = 5;
    mS = 0.8;
    mrk = 1;
    mr0 = 0.08;
    mdE = -0.07;
    mTrk = 14;
    
    % load all station data & find closest grid point
    Nst = 0;
    for i=1:N
        sds_i = load_station_data_for_wrf(S(i).name(1:5),Ts);
        sds_i.stid = S(i).name(1:5);
        [lati,loni,disti] = find_closest_point(sds_i.lat,sds_i.lon,glat,glon);
        if(disti < 3)
            Nst = Nst + 1;
            sds_i.glat = lati;
            sds_i.glon = loni;
            sds{Nst} = sds_i;
        else
            fprintf('Skipping station %s, too far from grid (%g km)\n', S(i).name(1:5), disti);
        end
    end
    sds = sds(1:Nst);
    
    % initialize all models
    ms = zeros(Nlon,Nlat,M);
    for i=1:Nlon
        for j=1:Nlat
            ms(i,j,1:3) = 0.5*(ew(i,j,1)+ed(i,j,1));
            ms(i,j,4) = mdE;
        end
    end        

    % diagnostic variables at station locations
    sts_per_time  = zeros(Nt,1);
    fm10_tgt      = nan * ones(Nt,Nst);
    fm10_model    = zeros(Nt,Nst,M);

    % for each timepoint in Ts
    for t=2:length(Ts)

        neg_fcast = zeros(3,1);
        
        % find stations that have fm10 observations at this time
        fm10o = nan * ones(Nst,1);
        fm10v = nan * ones(Nst,1);
        for i=1:Nst
            if(isfinite(sds{i}.fm10(t)))
                fm10o(i) = sds{i}.fm10(t);
                fm10_tgt(t,i) = fm10o(i);
                fm10v(i) = sds{i}.fm10_var(t);
            end
        end
        sts_per_time(t) = Nst;
        
        fprintf('time[%d]=%10.3f: sts=%d\n',t,Ts(t),Nst);

        % execute forecast procedure (for each filter) of the model to the
        % current time point using WRF environmental variables
        for i=1:Nlon
            for j=1:Nlat
                ed2 = 0.5 * (ed(i,j,t)+ed(i,j,t-1));
                ew2 = 0.5 * (ew(i,j,t)+ew(i,j,t-1));
                ri = rain(i,j,t);
                dt = (Ts(t)-Ts(t-1))*86400;
                mi = moisture_model_ext2(Tk,ed2,ew2,squeeze(ms(i,j,:)),ri,dt,1e10,mS,mrk,mr0,mTrk);
                neg_fcast = neg_fcast + (mi(1:3) < 0);
                mi(mi < 0) = 0;
                ms(i,j,:) = mi';
            end
        end
        
        if(any(neg_fcast>0))
            warning('Negative moistures [%d,%d,%d] found at in forecast at time [%d]=%g, corrected to zero',...
                neg_fcast(1),neg_fcast(2),neg_fcast(3),t,Ts(t));
        end
        
        % now store all diagnostic variables at station locations
        for s=1:Nst
            i = sds{s}.glon;
            j = sds{s}.glat;
            fm10_tgt(t,s) = fm10o(s);
            fm10_model(t,s,:) = ms(i,j,:);
        end
    end
            
    % dump the results to file
    save([out_dir,'/wrf_model_',out_dir,'.mat'], 'Ts','fm10_tgt','sds',...
         'fm10_model', 'rain', 'relh', 't2', 'ew', 'ed');
    
    % after all computations complete, show some pretty pictures`
    for s=1:Nst
        sd_s = sds{s};
        tsb = Ts - base_tday;
        i = sd_s.glon;
        j = sd_s.glat;
        rain_ij = squeeze(rain(i,j,:));
        ew_ij = squeeze(ew(i,j,:));
        
        f1 = figure();
        set(f1,'units','centimeters');
        set(f1,'papersize',[18,6])        
        set(f1,'paperposition',[0,0,18,6]);

        subplot(311);
        plot([tsb,tsb],[fm10_tgt(:,s),fm10_model(:,s,2)],'linewidth',1.2); 
        legend('tgt', 'model');
        title([sd_s.stid,': Observation tracking & forecasting']);
        ylabel('fm10 [-]');
        
        subplot(312);
        plot([tsb,tsb],[sd_s.rain,rain_ij],'linewidth',1.2);
        legend('st. rain','wrf rain');
        ylabel('Rain [mm/h]');

        subplot(313);
        plot([tsb,tsb],[sd_s.ews,ew_ij]);
        legend('station Ew', 'WRF Ew');
        ylabel('equilibrium [g/100g]');
        xlabel('Time [days from 1.1]');
        
        print('-dpng', [out_dir,'/',out_dir,'_',sd_s.stid,'_model_vs_time']);

        f2 = figure();
        set(f2,'units','centimeters');
        set(f2,'papersize',[12,12]);
        set(f2,'paperposition',[0,0,12,12]);        

        subplot(231);
        plot(fm10_tgt(:,s),fm10_model(:,s,2), 'o', 'MarkerSize',5);
        hold on;
        plot([0,max([fm10_tgt(:,s);fm10_model(:,s,2)])],[0,max([fm10_tgt(:,s);fm10_model(:,s,2)])],'k-');
        hold off;
        axis('equal');
        title([sd_s.stid,': Observation vs. Model fit']);
        xlabel('observation [-]');
        ylabel('model fit [-]');

        subplot(232);
        plot(sd_s.rain,rain_ij,'o','MarkerSize',5);
        hold on;
        plot([min([sd_s.rain;rain_ij]),max([sd_s.rain;rain_ij])], ...
             [min([sd_s.rain;rain_ij]),max([sd_s.rain;rain_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station rain vs. wrf rain']);
        axis('equal');
        xlabel('station rain [mm/h]');
        ylabel('WRF rain [mm/h]');

        subplot(233);
        t2_ij = squeeze(t2(i,j,:));
        plot(sd_s.t2,t2_ij,'o','MarkerSize',5);
        hold on;
        plot([min([sd_s.t2;t2_ij]),max([sd_s.t2;t2_ij])], ...
             [min([sd_s.t2;t2_ij]),max([sd_s.t2;t2_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station T2 vs. WRF T2']);
        axis('equal');
        xlabel('station t2 [K]');
        ylabel('WRF t2 [K]');

        subplot(234);
        relh_ij = squeeze(relh(i,j,:));
        plot(sd_s.relh,relh_ij,'o','MarkerSize',5);
        hold on;
        plot([0,max([sd_s.relh;relh_ij])],[0,max([sd_s.relh;relh_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station RELH vs. WRF computed RELH']);
        axis('equal');
        xlabel('station relh [-]');
        ylabel('wrf relh [-]');

        subplot(235);
        plot(sd_s.ews,ew_ij,'o','MarkerSize',5);
        hold on;
        plot([0,max([sd_s.ews;ew_ij])],[0,max([sd_s.ews;ew_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station E_w vs. WRF computed E_w']);
        axis('equal');
        xlabel('station E_w [-]');
        ylabel('wrf E_w [-]');

        print('-dpng', [out_dir,'/',out_dir,'_',sd_s.stid,'_model_scatters']);

        % report
        valids = find(isfinite(fm10_tgt(:,s)));
        Nt = length(valids);
        fprintf('Station: %s MODEL/MAPE: %g, MODEL/MSE: %g\n', ...
                sd_s.stid, norm(fm10_model(valids,s,2)-fm10_tgt(valids,s),1)/Nt, norm(fm10_model(valids,s,2)-fm10_tgt(valids,s),2)/Nt);
    end
        
end
