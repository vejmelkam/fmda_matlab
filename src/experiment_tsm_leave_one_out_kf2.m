

%
% This function conducts an experiment evaluating the ability of the TSM to
% transport observation information across space using the Fay-Herriott
% model.  The testing is done via a leave-one-out testing strategy.
%
%  synopsis: out = experiment_tsm_leave_one_out_kf(station_start,station_skip,filter_type)
%

function out = experiment_tsm_leave_one_out_kf2(station_start,station_skip,year,filter_type)

    % find all stations with data in 2012
    years = num2str(year);
    S = dir(['../data/*fm10*',years,'*']);
    N = length(S);

    base_tday = datenum(year,1,1,0,0,0);
    sds = cell(N,1);
    max_tday = 0;
    
    out = cell(N,1);

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
    
    % load all station data
    for i=1:N
        sds{i} = load_station_data(S(i).name(1:5),years);
        sds{i}.tdays = sds{i}.tdays - base_tday;
        max_tday = max(max_tday, sds{i}.tdays(end));
    end

    % Recompute the LOO procedure with each station left out in turn
    for z=station_start:station_skip:length(S)
        Ts = 120:1/24:min(270,max_tday);
        sts_per_time = zeros(numel(Ts),1);
        loo_ndx  = z;
        loo_tgt  = zeros(numel(Ts),1);
        loo_tsmc = zeros(numel(Ts),1);
        loo_tsm  = zeros(numel(Ts),1);
        loo_tm   = zeros(numel(Ts),1);
        loo_var  = zeros(numel(Ts),1);
        loo_varc = zeros(numel(Ts),1);
        betas    = zeros(numel(Ts),1);
        betasc   = zeros(numel(Ts),1);
        station_code = sds{z}.stid;
        fprintf('Leaving out station %d: %s\n', z,station_code);

        % re-initialize all models
        m_init = false(N,1);
        m_lastt = zeros(N,1);
        ms = zeros(N,M);
        P = zeros(N,M,M);
        sigma_fs = zeros(N,2*M+1,2*(2*M+1)+1);
        sqrtP = zeros(N,M,2*(2*M+1)+1);
        for i=1:N
            ndxs = find_valid_obs(120,sds(i));
            P(i,:,:) = 0.01 * eye(M);
            if(isfinite(ndxs))
                ms(i,1:3) = 0.5*(sds{i}.ew(ndxs(1))+sds{i}.ed(ndxs(1)));
                ms(i,4) = -0.04;
                m_init(i) = true;
                m_lastt(i) = 120;
            end
        end        
        fprintf('Initialized %d models out of %d\n', sum(m_init), N);

        % for each (hourly) timepoint in Ts
        iter_ndx = 1;
        conv_failures = 0;
        for t=2:length(Ts)

            % find stations that have (rel_humidity,air_temp,accum_precip)
            % observations valid for this time.  Additionally, read off their
            % fm10 observations (which still may be nans)
            [ndxs,fm10o,fm10v] = find_valid_obs(Ts(t),sds);
            st_ndx = find(isfinite(ndxs));
            Nst = length(st_ndx);
            sts_per_time(t) = Nst;
            times = ndxs(st_ndx);
            st_map = false(Nst,1);

            if(rem(t,100)==0)
                fprintf('%d ', t);
            end

            % if the model has an observation valid now, integrate the model
            % to this time point (follow UKF forecast procedure)
            % mark model as updated to this time point
            for si=1:Nst
                i = st_ndx(si);
                if(m_init(i))
                    ed = sds{i}.ed(ndxs(i));
                    ew = sds{i}.ew(ndxs(i));
                    ri = sds{i}.rain(ndxs(i));
                    dt = (Ts(t)-m_lastt(i))*86400;
                    if(dt<=0)
                        error('dt must be positive')
                    end
                    Pi = squeeze(P(i,:,:));
                    if(filter_type==2)
                        f = @(x,w) moisture_model_ext2(Tk,ed,ew,x,ri,dt,1e10,0.6,2,0.08,7) + w;
                        [mi,sqrtPi,sigma_f] = ukf_forecast_general(ms(i,:)',f,Pi,Qphr*dt/3600,1,kappa);
                        sigma_fs(i,:,:) = sigma_f;
                        sqrtP(i,:,:) = sqrtPi;
                        Pi = sqrtPi*sqrtPi';
                    elseif(filter_type==1)
                        [mi,Pi] = ekf_forecast2(Tk,ed,ew,ms(i,:)',ri,dt,1e10,Pi,Qphr,0.6,2,0.08,7);
                    else
                        error('Invalid filter type');
                    end
                    if(any(isnan(Pi)))
                        warning('nans found in forecast covariance');
                        Pi
                    end
                    if(any(eig(Pi) < 0))
                        Pi
                        warning('negative eigenvalues in forecast covariance');
                    end
                    ms(i,:) = mi';
                    P(i,:,:) = Pi;
                    m_lastt(i) = Ts(t);
                    st_map(si) = true;
                else
                    ms(i,1:3) = 0.5*(sds{i}.ew(ndxs(i))+sds{i}.ed(ndxs(i)));
                    ms(i,4) = -0.04;
                    m_init(i) = true;
                    m_lastt(i) = Ts(t);
                    st_map(si) = false;
                end
            end
            
            % check for negative moistures
            if(any(ms(:,1:3) < 0))
                warning('Negative moistures found!');
                ms
            end

            % construct regressors, observations and variances: only use models
            % that have been integrated to this time point (marked in st_map).
            X = zeros(Nst,4);
            Z = zeros(Nst,1); 
            G = zeros(Nst,1);
            Nobs = 0;
            for o=1:Nst
                so = st_ndx(o);
                if(st_map(o) && isfinite(fm10o(so)))
                    std = sds{so};
                    stt = times(o);
                    Nobs = Nobs + 1;
                    X(Nobs,:) = [ms(so,2),1,std.elevation/2000,std.rain(stt)];
                    Z(Nobs) = fm10o(so);
                    G(Nobs) = fm10v(so);
                end
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
                loo_in_st = find(loo_ndx == st_ndx);
                if(~isempty(loo_in_st))
                    loo_tm(iter_ndx) = Ts(t);
                    loo_sd = sds{loo_ndx};
                    stt = times(loo_in_st);
                    loo_tgt(iter_ndx) = loo_sd.fm10(stt);
                    y_loo = [ms(loo_ndx,2),1,loo_sd.elevation/2000,loo_sd.rain(stt)]';
                    Xe = [X;y_loo(1:size(X,2))'];
                else
                    Xe = X;
                end
                [beta,sigma2] = estimate_tsm_parameters(X,Z,G);
                [betac,sigma2c] = estimate_tsm_parameters_constr(X,Z,G,Xe,0.6);
                if(isnan(sigma2c))
                    sigma2c = sigma2;
                    betac = beta;
                    conv_failures = conv_failures + 1;
                end
                betas(iter_ndx) = beta(1);
                betasc(iter_ndx) = betac(1);

                % compute X'Sigma^{-1}*X
                XSX = (X'*diag(1./(sigma2+G))*X);
                XSXc = (X'*diag(1./(sigma2c+G))*X);

                if(rcond(XSXc) < 1e-16)
                    warning('XSXc is badly conditioned!');
                    [X,Z,G]
                end
                % do the following (LOO testing) only if the left-out stations
                % actually has at least valid weather observations here
                % fm10 can still be unavailable, in which case loo_tgt(iter_ndx)
                % becomes a nan
                loo_in_st = find(loo_ndx == st_ndx);
                loo_tm(iter_ndx) = Ts(t);
                if(~isempty(loo_in_st))
                    loo_sd = sds{loo_ndx};
                    loo_tgt(iter_ndx) = loo_sd.fm10(times(loo_in_st));
                    stt = times(loo_in_st);
                    y_loo = [ms(loo_ndx,2),1,loo_sd.elevation/2000,loo_sd.rain(stt)]';
                    y_loo = y_loo(1:size(X,2));
                    loo_tsm(iter_ndx) = y_loo' * beta;
                    loo_var(iter_ndx) = sigma2 + y_loo' * (XSX\y_loo);
                    loo_tsmc(iter_ndx) = y_loo' * betac;
                    loo_varc(iter_ndx) = sigma2c + y_loo' * (XSXc\y_loo);
                else
                    loo_tgt(iter_ndx) = nan;
                end

                % find regressors for all included stations
                for o=1:Nst
                    if(st_map(o))
                        so = st_ndx(o);
                        std = sds{so};
                        stt = times(o);
                        xr = [ms(so,2),1,std.elevation/2000,std.rain(stt)]';
                        xr = xr(1:size(X,2));

                        % remove the part explained by the forecast
%                        xr(1) = 0;
                        xpseudo = xr' * betac;
                        varpseudo = sigma2c + xr' * (XSXc\xr);
                        if(filter_type==2)
                            sqrtPo = squeeze(sqrtP(so,:,:));
                            sigma_f = squeeze(sigma_fs(so,:,:));
                            [ms(so,:),P(so,:,:)] = ukf_update(ms(so,:)',sqrtPo,sigma_f,[0,1,0,0,0],xpseudo,varpseudo,kappa);
                        elseif(filter_type==1)
                            Po = squeeze(P(so,:,:));
                            [ms(so,:),P(so,:,:)] = ekf_update(ms(so,:)',Po,[0,1,0,0,0],xpseudo,varpseudo);
                        else
                            error('Invalid filter type');
                        end
                        
                        if(any(isnan(P(so,:,:))))
                            squeeze(P(so,:,:))
                            error('nans found in updated covariance');
                        end
                        if(any(eig(squeeze(P(so,:,:))) < 0))
                            squeeze(P(so,:,:))
                            error('negative eigenvalues found in updated covariance');
                        end
                    end
                end
                
                % check for negative moistures
                if(any(ms(:,1:3) < 0))
                    warning('negative moisture post-assimilation');
                    ms
                end
                
                % we have had an assimilation
                iter_ndx = iter_ndx + 1;
            else
                fprintf('No observations at time %g\n', Ts(t));
            end

        end

        % we censor loo_tsm to be between 0 and 0.6
        loo_tsm(loo_tsm < 0) = 0;
        loo_tsm(loo_tsm > 0.6) = 0.6;

        loo_tgt = loo_tgt(1:iter_ndx-1);
        loo_tsm = loo_tsm(1:iter_ndx-1);
        loo_tsmc = loo_tsmc(1:iter_ndx-1);
        loo_tm = loo_tm(1:iter_ndx-1);
        loo_var = loo_var(1:iter_ndx-1);
        loo_varc = loo_varc(1:iter_ndx-1);
        betas = betas(1:iter_ndx-1);
        betasc = betasc(1:iter_ndx-1);

        if(iter_ndx > 5)

            f1 = figure();
            set(f1,'units','centimeters');
            set(f1,'papersize',[18,6])        
            set(f1,'paperposition',[0,0,18,6]);

            subplot(311);
            plot([loo_tm,loo_tm,loo_tm],[loo_tgt,loo_tsm,loo_tsmc],'linewidth',1.2);
            legend('tgt', 'tsm', 'tsmc');
            title([sds{z}.stid,': Comparison between TSM fit (LS and CLS) and observation']);
            xlabel('Time [days from 1.1]');
            ylabel('fm10 [-]');

            subplot(312);
            plot([loo_tm,loo_tm],[abs(loo_tgt-loo_tsm),sqrt(loo_var)]);
            title([sds{z}.stid,': TSM-obs error vs. LS/TSM variance estimate']);
            legend('ls: abs err', 'ls: stdev');
            xlabel('Time [days from 1.1]');
            ylabel('fm10 [-]');

            subplot(313);
            plot([loo_tm,loo_tm],[abs(loo_tgt-loo_tsmc),sqrt(loo_varc)]);
            title([sds{z}.stid,': TSM-obs error vs. CLS/TSM variance estimate']);
            legend('cls: abs err', 'cls: stdev');
            xlabel('Time [days from 1.1]');
            ylabel('fm10 [-]');

            print('-dpng', [station_code,'_',years,'_tsm2_vs_time_kf']);

            f2 = figure();
            set(f2,'units','centimeters');
            set(f2,'papersize',[12,12]);
            set(f2,'paperposition',[0,0,12,12]);        

            subplot(221);
            plot(loo_tgt,loo_tsm, 'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
            hold on;
            plot([0,0.5],[0,0.5],'b-');
            hold off;
            axis('square');
            title([sds{z}.stid,': Observation vs. LS TSM fit']);
            xlabel('observation [-]');
            ylabel('TSM fit [-]');

            subplot(222);
            plot(abs(loo_tgt-loo_tsm),sqrt(loo_var),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
            title([sds{z}.stid,': Abs. error vs. LS TSM stdev estimate']);
            axis('square');
            xlabel('abs. error [-]');
            ylabel('stdev [-]');

            subplot(223);
            plot(loo_tgt,loo_tsmc, 'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
            hold on;
            plot([0,0.5],[0,0.5],'b-');
            hold off;
            axis('square');
            title([sds{z}.stid,': Observation vs. CLS TSM fit']);
            xlabel('observation [-]');
            ylabel('TSM fit [-]');

            subplot(224);
            plot(abs(loo_tgt-loo_tsmc),sqrt(loo_varc),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5);
            title([sds{z}.stid,': Abs. error vs. CLS TSM stdev estimate']);
            axis('square');
            xlabel('abs. error [-]');
            ylabel('stdev [-]');

            print('-dpng', [station_code,'_',years,'_tsm2_scatters_kf']);

        end

        % report
        valids = find(isfinite(loo_tgt));
        Nt = length(valids);
        fprintf('\n\n*** report ***\n');
        fprintf('Station: %s LS/MAPE: %g  LS/MSE: %g  CLS/MAPE: %g   CLS/MSE: %g  convergence failures: %d\n', ...
                station_code, norm(loo_tgt(valids)-loo_tsm(valids),1)/Nt, norm(loo_tgt(valids)-loo_tsm(valids),2)/Nt, ...
                norm(loo_tgt(valids)-loo_tsmc(valids),1)/Nt, norm(loo_tgt(valids)-loo_tsmc(valids),2)/Nt, conv_failures);
        
        % dump the results to file
        sdsz = sds{z};
        save([station_code,'_',years,'_tsm2_loo.mat'],'loo_tm','loo_tsm','loo_tsmc','loo_tgt','betas', 'betasc','sdsz','loo_var','loo_varc');

    end
end
