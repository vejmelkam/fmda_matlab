

%
% This script conducts an experiment evaluating the ability of the TSM to
% transport observation information across space using the Fay-Herriott
% model.  This is done via a leave-one-out testing strategy.
%

S = dir('../data/*fm10*2012*');
N = length(S);

year = '2012';
base_tday = datenum(2012,1,1,0,0,0);
sds = cell(N,1);
max_tday = 0;

Tk = [1,10,100]';
M = 5;
k = 3;

Qphr = zeros(M);
Qphr(1:k, 1:k) = diag([0.0005,0.00005,0.00001]);
Qphr(k+1,k+1) = 0.0001;
Qphr(k+2,k+2) = 0.0001;

for i=1:N
    sds{i} = load_station_data(S(i).name(1:5),year);
    sds{i}.tdays = sds{i}.tdays - base_tday;
    max_tday = max(max_tday, sds{i}.tdays(end));
end


for z=2:length(S)
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
    fprintf('Leaving out station: %s\n', station_code);
    
    % re-initialize all models
    m_init = false(N,1);
    m_lastt = zeros(N,1);
    ms = zeros(N,M);
    P = zeros(N,M,M);
    sqrtP = zeros(N,M,2*M+1);
    for i=1:N
        ndxs = find_valid_obs(120,sds(i));
        P(i,:,:) = 0.01 * eye(M);
        if(isfinite(ndxs))
            ms(i,1:3) = 0.5*(sds{i}.ew(ndxs(1))+sds{i}.ed(ndxs(1)));
            m_init(i) = true;
            m_lastt(i) = 120;
        end
    end        
    fprintf('Initialized %d models out of %d\n', sum(m_init), N);
    
    iter_ndx = 1;
    for t=1:length(Ts)

        [ndxs,fm10o,fm10v] = find_valid_obs(Ts(t),sds);
        st_ndx = find(isfinite(ndxs));
        Nst = length(st_ndx);
        sts_per_time(t) = Nst;
        times = ndxs(st_ndx);
        st_map = false(Nst,1);
        
        if(rem(t,250)==0)
            fprintf('t=%d\n', t);
        end
        
        % run each model to Ts if possible
        for si=1:length(st_ndx)
            i = st_ndx(si);
            if(m_init(i))
                ed = sds{i}.ed(ndxs(i));
                ew = sds{i}.ew(ndxs(i));
                ri = sds{i}.rain(ndxs(i));
                dt = (Ts(t)-m_lastt(i))*86400;
                Pi = squeeze(P(i,:,:));
                if(any(isnan(Pi(:))))
                    Pi
                end
                [mi,sqrtPi,Pi] = ukf_forecast(Tk,ed,ew,ms(i,:)',ri,dt,1e10,Pi,Qphr);
                if(any(isnan(Pi(:))))
                    Pi
                end
                if(any(isnan(sqrtP(:))))
                    sqrtP
                end
                ms(i,:) = mi';
                sqrtP(i,:,:) = sqrtPi;
                P(i,:,:) = Pi;
                m_lastt(i) = Ts(t);
                st_map(si) = true;
            else
                ms(i,1:3) = 0.5*(sds{i}.ew(ndxs(i))+sds{i}.ed(ndxs(i)));
                m_init(i) = true;
                m_lastt(i) = Ts(t);
                st_map(si) = true;
            end
        end

        % construct regressors, observations and variances
%        fprintf('Time %g - have %d valid obs building %d regressors\n', Ts(t),length(st_ndx),sum(reg_map));
        X = zeros(Nst,4);
        Z = zeros(Nst,1); 
        G = zeros(Nst,1);
        ndx = 1;
        for o=1:Nst
            if(st_map(o) && isfinite(fm10o(st_ndx(o))))
                std = sds{st_ndx(o)};
                stt = times(o);
                X(ndx,:) = [ms(st_ndx(o),2),1,std.elevation,std.rain(stt)];
                Z(ndx) = fm10o(st_ndx(o));
                G(ndx) = fm10v(st_ndx(o));
                ndx = ndx + 1;
            end
        end
        if(all(X(:,4)==0))
            X = X(:,1:3);
        end
        X = X(1:ndx-1,:);
        Z = Z(1:ndx-1,:);
        G = G(1:ndx-1,:);
        
        if(size(X,2) > size(X,1))
            X = X(:,1:size(X,1));
        end            
        
        [beta,sigma2] = estimate_tsm_parameters(X,Z,G);
        [betac,sigma2c] = estimate_tsm_parameters_constr(X,Z,G,X);
        if(isnan(sigma2c))
            sigma2c = sigma2;
            betac = beta;
        end
        betas(iter_ndx) = beta(1);
        betasc(iter_ndx) = betac(1);
        
        % find regressor vector for target station and make prediction
        XSX = (X'*diag(1./(sigma2+G))*X);
        XSXc = (X'*diag(1./(sigma2c+G))*X);
        % find our left-out station in the array
        loo_in_st = find(loo_ndx == st_ndx);
        loo_tm(iter_ndx) = Ts(t);
        if(~isempty(loo_in_st))
            loo_sd = sds{loo_ndx};
            loo_tgt(iter_ndx) = loo_sd.fm10(times(loo_in_st));
            stt = times(loo_in_st);
            y_loo = [ms(loo_ndx,2),1,loo_sd.elevation,loo_sd.rain(stt)]';
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
                std = sds{st_ndx(o)};
                stt = times(o);
                xr = [ms(st_ndx(o),2),1,std.elevation,std.rain(stt)]';
                xr = xr(1:size(X,2));
                % pseudo-observations are censored
                xpseudo = min(max(xr' * betac,0),2.5);
                varpseudo = sigma2c + xr' * (XSXc\xr) + 1e-4;
                Po = squeeze(P(o,:,:));
                sqrtPo = squeeze(sqrtP(o,:,:));
                [ms(o,:),P(o,:,:)] = ukf_update(ms(o,:)',sqrtPo,Po,[0,1,0,0,0],xpseudo,varpseudo);
                if(any(isnan(P(o,:,:))))
                    P(o,:,:)
                end
            end
        end

        iter_ndx = iter_ndx + 1;

    end
    
    % we censor loo_tsm to be between 0 and 2.5
    loo_tsm(loo_tsm < 0) = 0;
    loo_tsm(loo_tsm > 2.5) = 2.5;

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
        
        print('-dpng', [station_code,'_tsm_vs_time_kf']);

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
        
        print('-dpng', [station_code,'_tsm_scatters_kf']);
        
        close all;

    end
    
    % report
    valids = find(isfinite(loo_tgt));
    Nt = length(valids);
    fprintf('Station: %s LS/MAPE: %g  LS/MSE: %g  CLS/MAPE: %g   CLS/MSE: %g\n', ...
            station_code, norm(loo_tgt(valids)-loo_tsm(valids),1)/Nt, norm(loo_tgt(valids)-loo_tsm(valids),2)/Nt, ...
            norm(loo_tgt(valids)-loo_tsmc(valids),1)/Nt, norm(loo_tgt(valids)-loo_tsmc(valids),2)/Nt);
        
end

