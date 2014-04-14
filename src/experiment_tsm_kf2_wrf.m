

%
% This function conducts an experiment evaluating the ability of the TSM to
% transport observation information across space using the Fay-Herriott
% model.  The testing is done via a leave-one-out testing strategy.
%
%  synopsis: out = experiment_tsm_leave_one_out_kf(station_start,station_skip,filter_type)
%

function experiment_tsm_kf2_wrf(wrf_file,filter_type,stop_da_at)

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
    hgt = ncread(wrf_file,'HGT');
    hgt = hgt(:,:,1);
    rain = cat(3,zeros(size(glon)),diff(ncread(wrf_file,'RAINC'),1,3) + diff(ncread(wrf_file,'RAINNC'),1,3));
    [Nlon,Nlat] = size(glon);
    [ed,ew,relh] = equilibrium_moisture(psfc,q2,t2);
    relh = relh / 100;
    
    assert(all(ed(:) > 0));
    assert(all(ew(:) > 0));

    relh_neg = sum(relh(:) < 0);
    relh_pos = sum(relh(:) > 1);
    if(relh_neg > 0)
        warning('fixing WRF relative humidity, found %d negative values', relh_neg);
        relh(relh < 0) = 0;
    end
    if(relh_pos > 0)
        warning('fixing WRF relative humidity, found %d values over 1', relh_pos);
        relh(relh > 1) = 1;
    end
    
    % remove any rain that is too light to affect the moisture model
    rain(rain < 0.08) = 0;
    
    % parse the dates
    Ts = zeros(Nt,1);
    for i=1:Nt
        ti = wrf_t(i,:);
        Ts(i) = datenum(str2double(ti(1:4)),str2double(ti(6:7)),str2double(ti(9:10)),...
                        str2double(ti(12:13)),str2double(ti(15:16)),str2double(ti(18:19)));
        % convert rain to mm/h units
        if(i>1)
            rain(:,:,i) = rain(:,:,i) / (Ts(i) - Ts(i-1)) * 3600;
        end
    end
    
    % load data from correct year
    year = str2double(ti(1:4));
    years = num2str(year);
    S = dir(['../data/*fm10*',years,'*']);
    N = length(S);

    sds = cell(N,1);
    base_tday = datenum(year,1,1,0,0,0);
    
    % output dir prefix
    if(filter_type==2), flt_suffix = 'ukf'; else flt_suffix='ekf'; end
    ti = wrf_t(1,:);
    out_dir = sprintf('wrf_%d%d%d_%s', str2double(ti(1:4)), str2double(ti(6:7)), str2double(ti(9:10)), flt_suffix);
    mkdir(out_dir);
    
    % constant model parameters
    Tk = [1,10,100]';
    M = 5;
    k = 3;

    % process noise matrix
    Qphr = zeros(M);
    Qphr(1:k, 1:k) = diag([0.0005,0.00005,0.00001]);
    Qphr(k+1,k+1) = 0.0001;
    Qphr(k+2,k+2) = 0.0001;
    
    kappa = 0;
    
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
    P = zeros(Nlon,Nlat,M,M);
    sigma_fs = zeros(Nlon,Nlat,2*M+1,2*(2*M+1)+1);
    sqrtP = zeros(Nlon,Nlat,M,2*(2*M+1)+1);
    for i=1:Nlon
        for j=1:Nlat
            P(i,j,:,:) = 0.01 * eye(M);
            ms(i,j,1:3) = 0.5*(ew(i,j,1)+ed(i,j,1));
            ms(i,j,4) = -0.04;
        end
    end        

    % diagnostic variables at station locations
    sts_per_time  = zeros(Nt,1);
    fm10_tgt      = nan * ones(Nt,Nst);
    fm10_tsmc     = zeros(Nt,Nst);
    fm10_tsm      = zeros(Nt,Nst);
    fm10_tsm_var  = zeros(Nt,Nst);
    fm10_tsmc_var = zeros(Nt,Nst);
    fm10_model    = zeros(Nt,Nst,M);
    betas         = zeros(Nt,1);
    betasc        = zeros(Nt,1);

    % for each timepoint in Ts
    conv_failures = 0;
    for t=2:length(Ts)

        neg_tsmc = 0;
        neg_tsm = 0;
        neg_fcast = zeros(3,1);
        
        % find stations that have fm10 observations at this time
        valid_now = false(Nst,1);
        fm10o = nan * ones(Nst,1);
        fm10v = nan * ones(Nst,1);
        for i=1:Nst
            if(isfinite(sds{i}.fm10(t)))
                valid_now(i) = true;
                fm10o(i) = sds{i}.fm10(t);
                fm10_tgt(t,i) = fm10o(i);
                fm10v(i) = sds{i}.fm10_var(t);
            end
        end
        st_ndx = find(valid_now)';
        Nvalidnow = length(st_ndx);
        sts_per_time(t) = Nst;
        
        fprintf('time[%d]=%10.3f: sts=%d\n',t,Ts(t),Nst);

        % execute forecast procedure (for each filter) of the model to the
        % current time point using WRF environmental variables
        for i=1:Nlon
            for j=1:Nlat
                ed2 = 0.5 * (ed(i,j,t)+ed(i,j,t-1));
                ew2 = 0.5 * (ew(i,j,t)+ew(i,j,t-1));
                assert(ed2 > 0);
                assert(ew2 > 0);
                ri = rain(i,j,t);
                dt = (Ts(t)-Ts(t-1))*86400;
                Pi = squeeze(P(i,j,:,:));
                if(filter_type==2)
                    f = @(x,w) moisture_model_ext2(Tk,ed2,ew2,x,ri,dt,1e10,0.6,2,0.08,7) + w;
                    [mi,sqrtPi,sigma_f] = ukf_forecast_general(squeeze(ms(i,j,:)),f,Pi,Qphr*dt/3600,1,kappa);
                    sigma_fs(i,j,:,:) = sigma_f;
                    sqrtP(i,j,:,:) = sqrtPi;
                    Pi = sqrtPi*sqrtPi';
                elseif(filter_type==1)
                    [mi,Pi] = ekf_forecast2(Tk,ed2,ew2,squeeze(ms(i,j,:)),ri,dt,1e10,Pi,Qphr,0.6,2,0.08,7);
                else
                    error('Invalid filter type');
                end
                if(any(isnan(Pi)))
                    Pi
                    error('nans found in forecast covariance at (%d,%d)',i,j);
                end
                if(any(eig(Pi) < 0))
                    Pi
                    error('negative eigenvalues in forecast covariance at (%d,%d)',i,j);
                end
                neg_fcast = neg_fcast + (mi(1:3) < 0);
                mi(mi < 0) = 0;
                ms(i,j,:) = mi';
                P(i,j,:,:) = Pi;
            end
        end
        
        if(any(neg_fcast>0))
            warning('Negative moistures [%d,%d,%d] found at in forecast at time [%d]=%g, corrected to zero',...
                neg_fcast(1),neg_fcast(2),neg_fcast(3),t,Ts(t));
        end
        
        % construct regressors, observations and variances: only use models
        % that have been integrated to this time point (marked in st_map).
        X = zeros(Nvalidnow,4);
        Z = zeros(Nvalidnow,1); 
        G = zeros(Nvalidnow,1);
        Nobs = 0;
        for o=st_ndx
            std = sds{o};
            i = std.glon;
            j = std.glat;
            Nobs = Nobs + 1;
            X(Nobs,:) = [ms(i,j,2),1,hgt(i,j)/2000,rain(i,j,t)];
            Z(Nobs) = fm10o(o);
            G(Nobs) = fm10v(o);
            fm10_tgt(t,o) = fm10o(o);
        end

        % the fourth regressor is rain intensity - if that is zero
        % everywhere, we remove the regressor
        if(all(X(:,4)==0))
            X = X(:,1:3);
        end

        % only use the stations that we filled out (ndx-1)
        X = X(1:Nobs,:);
        Z = Z(1:Nobs,:);
        G = G(1:Nobs,:);

        % if there are less observations than regressors, reduce the number of
        % regressors to the number of observations
        if(size(X,2) > size(X,1))
            X = X(:,1:size(X,1));
        end            

        % run both methods - the constraint tsm and the unconstrained tsm
        % if the constrained method does not converge, then replace its
        % result by the unconstrained solver (as it would in an operational
        % system)
        if(Nobs > 0)
            Xe = zeros(size(X,1)+400,4);
            for z=size(X,1)+1:size(Xe,1)
                i = randi(Nlon);
                j = randi(Nlat);
                Xe(z,:) = [ms(i,j,2),1,hgt(i,j)/2000,rain(i,j,t)];
            end
            Xe = Xe(:,1:size(X,2));
            [beta,sigma2] = estimate_tsm_parameters(X,Z,G);
           [betac,sigma2c] = estimate_tsm_parameters_constr(X,Z,G,Xe,0.6);
           if(isnan(sigma2c))
               sigma2c = sigma2;
               betac = beta;
               conv_failures = conv_failures + 1;
           end
           betas(t) = beta(1);
           betasc(t) = betac(1);

            % compute X'Sigma^{-1}*X
            XSX = (X'*diag(1./(sigma2+G))*X);
           XSXc = (X'*diag(1./(sigma2c+G))*X);

            if(rcond(XSX) < 1e-16)
                warning('XSX is badly conditioned!');
                [X,Z,G]
            end

            % run update step on all stations
            if(t <= stop_da_at)
                for i=1:Nlon
                    for j=1:Nlat
                        xr = [ms(i,j,2),1,hgt(i,j)/2000,rain(i,j,t)]';
                        xr = xr(1:size(X,2));
                        if(xr' * beta < 0), neg_tsm = neg_tsm + 1; end
                       if(xr' * betac < 0), neg_tsmc = neg_tsmc + 1; end
                        
                        % run 
                        xpseudo = min(max(xr' * beta,0),0.6);
                        varpseudo = sigma2 + xr' * (XSX\xr);
                        
                        if(filter_type==2)
                            sqrtPo = squeeze(sqrtP(i,j,:,:));
                            sigma_f = squeeze(sigma_fs(i,j,:,:));
                            [ms(i,j,:),P(i,j,:,:)] = ukf_update(squeeze(ms(i,j,:)),sqrtPo,sigma_f,[0,1,0,0,0],xpseudo,varpseudo,kappa);
                        elseif(filter_type==1)
                            Po = squeeze(P(i,j,:,:));
                            [ms(i,j,:),P(i,j,:,:)] = ekf_update(squeeze(ms(i,j,:)),Po,[0,1,0,0,0],xpseudo,varpseudo);
                        else
                            error('Invalid filter type');
                        end

                    if(any(isnan(P(i,j,:,:))))
                        squeeze(P(i,j,:,:))
                        error('nans found in updated covariance at (%d,%d)',i,j);
                    end
                    if(any(eig(squeeze(P(i,j,:,:))) < 0))
                        squeeze(P(i,j,:,:))
                        error('negative eigenvalues found in updated covariance at (%d,%d)',i,j);
                    end
                    end
                    % check for negative moistures
                    if(any(ms(i,j,1:3) < 0))
                        ms
                        warning('negative moisture post-assimilation at (%d,%d) at time [%d] %g',i,j,t,Ts(t));
                    end
                end
            end
            
            % now store all diagnostic variables at station locations
            for s=1:Nst
                i = sds{s}.glon;
                j = sds{s}.glat;
                xr = [ms(i,j,2),1,hgt(i,j)/2000,rain(i,j,t)]';
                xr = xr(1:size(X,2));

                fm10_tsmc(t,s) = xr' * betac;
                fm10_tsmc_var(t,s) = sigma2c + xr' * (XSXc\xr);

                fm10_tsm(t,s) = min(max(xr' * beta,0),0.6);
                fm10_tsm_var(t,s) = sigma2 + xr' * (XSX\xr);

                fm10_model(t,s,:) = ms(i,j,:);
            end
        end
        
        if((neg_tsmc > 0) || (neg_tsm > 0))
            warning('found negative tsm [%d] tsmc [%d] at time [%d] %g',neg_tsm,neg_tsmc,t,Ts(t));
        end
    end
            
    % dump the results to file
    save([out_dir,'/wrf_tsm2_',out_dir,'.mat'], ...
         'Ts','fm10_tgt','fm10_tsm','fm10_tsmc','fm10_tsm_var', ...
         'fm10_tsmc_var', 'betas', 'betasc','sds','fm10_model', ...
         'rain', 'relh', 't2', 'ew', 'ed', 'relh_neg','relh_pos');
    
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
        plot([tsb,tsb,tsb,tsb],[fm10_tgt(:,s),fm10_tsm(:,s),fm10_tsmc(:,s),fm10_model(:,s,2)],'linewidth',1.2); 
        hold on;
        plot([tsb(stop_da_at),tsb(stop_da_at)],[0,max(fm10_tgt(:,s))],'k-','linewidth',2);
        hold off;
        legend('tgt', 'tsm', 'tsmc', 'model');
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
        
        print('-dpng', [out_dir,'/',out_dir,'_',sd_s.stid,'_tsm2_vs_time_',flt_suffix]);

        f2 = figure();
        set(f2,'units','centimeters');
        set(f2,'papersize',[12,12]);
        set(f2,'paperposition',[0,0,12,12]);        

        subplot(331);
        plot(fm10_tgt(:,s),fm10_tsm(:,s), 'o', 'MarkerSize',5);
        hold on;
        plot([0,max([fm10_tgt(:,s);fm10_tsm(:,s)])],[0,max([fm10_tgt(:,s);fm10_tsm(:,s)])],'b-');
        hold off;
        axis('equal');
        title([sd_s.stid,': Observation vs. LS TSM fit']);
        xlabel('observation [-]');
        ylabel('TSM fit [-]');

        subplot(332);
        plot(fm10_tgt(:,s),fm10_tsmc(:,s), 'o', 'MarkerSize',5);
        hold on;
        plot([0,max([fm10_tgt(:,s);fm10_tsmc(:,s)])],[0,max([fm10_tgt(:,s);fm10_tsmc(:,s)])],'k-');
        hold off;
        axis('equal');
        title([sd_s.stid,': Observation vs. CLS TSM fit']);
        xlabel('observation [-]');
        ylabel('TSMc fit [-]');

        subplot(333);
        plot(fm10_tgt(:,s),fm10_model(:,s,2), 'o', 'MarkerSize',5);
        hold on;
        plot([0,max([fm10_tgt(:,s);fm10_model(:,s,2)])],[0,max([fm10_tgt(:,s);fm10_model(:,s,2)])],'k-');
        hold off;
        axis('equal');
        title([sd_s.stid,': Observation vs. Model fit']);
        xlabel('observation [-]');
        ylabel('model fit [-]');

        subplot(334);
        plot(abs(fm10_tgt(:,s)-fm10_tsm(:,s)),sqrt(fm10_tsm_var(:,s)),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
        title([sd_s.stid,': Abs. error vs. LS TSM stdev estimate']);
        axis('equal');
        xlabel('abs. error [-]');
        ylabel('stdev [-]');

        subplot(335);
        plot(abs(fm10_tgt(:,s)-fm10_tsmc(:,s)),sqrt(fm10_tsmc_var(:,s)),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
        title([sd_s.stid,': Abs. error vs. CLS TSM stdev estimate']);
        axis('equal');
        xlabel('abs. error [-]');
        ylabel('stdev [-]');

        subplot(336);
        plot(sd_s.rain,rain_ij,'o','MarkerSize',5);
        hold on;
        plot([min([sd_s.rain;rain_ij]),max([sd_s.rain;rain_ij])], ...
             [min([sd_s.rain;rain_ij]),max([sd_s.rain;rain_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station rain vs. wrf rain']);
        axis('equal');
        xlabel('station rain [mm/h]');
        ylabel('WRF rain [mm/h]');

        subplot(337);
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

        subplot(338);
        relh_ij = squeeze(relh(i,j,:));
        plot(sd_s.relh,relh_ij,'o','MarkerSize',5);
        hold on;
        plot([0,max([sd_s.relh;relh_ij])],[0,max([sd_s.relh;relh_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station RELH vs. WRF computed RELH']);
        axis('equal');
        xlabel('station relh [-]');
        ylabel('wrf relh [-]');

        subplot(339);
        plot(sd_s.ews,ew_ij,'o','MarkerSize',5);
        hold on;
        plot([0,max([sd_s.ews;ew_ij])],[0,max([sd_s.ews;ew_ij])],'k-');
        hold off;
        title([sd_s.stid,': Station E_w vs. WRF computed E_w']);
        axis('equal');
        xlabel('station E_w [-]');
        ylabel('wrf E_w [-]');

        print('-dpng', [out_dir,'/',out_dir,'_',sd_s.stid,'_tsm2_scatters_',flt_suffix]);

        % report
        valids = find(isfinite(fm10_tgt(:,s)));
        Nt = length(valids);
        fprintf('Station: %s MODEL/MAPE: %g, MODEL/MSE: %g LS/MAPE: %g  LS/MSE: %g  CLS/MAPE: %g   CLS/MSE: %g\n', ...
                sd_s.stid, norm(fm10_model(valids,s,2)-fm10_tgt(valids,s),1)/Nt, norm(fm10_model(valids,s,2)-fm10_tgt(valids,s),2)/Nt, ...
                norm(fm10_tgt(valids,s)-fm10_tsm(valids,s),1)/Nt, norm(fm10_tgt(valids,s)-fm10_tsm(valids,s),2)/Nt, ...
                norm(fm10_tgt(valids,s)-fm10_tsmc(valids,s),1)/Nt, norm(fm10_tgt(valids,s)-fm10_tsmc(valids,s),2)/Nt);
    end
        
end
