

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

for i=1:N
    sds{i} = load_station_data(S(i).name(1:5),year);
    sds{i}.tdays = sds{i}.tdays - base_tday;
    max_tday = max(max_tday, sds{i}.tdays(end));
end


for z=1:length(S)
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
    
    iter_ndx = 1;
    for t=1:length(Ts)

        [ndxs,fm10o] = find_valid_obs(Ts(t),sds);
        ndxs(isnan(fm10o)) = nan;
        st_ndx = find(isfinite(ndxs));
        Nst = length(st_ndx);
        sts_per_time(t) = Nst;
        times = ndxs(st_ndx);

        if((length(st_ndx) > 30) && any(st_ndx == loo_ndx))
            % find our left-out station in the array
            loo_in_st = find(loo_ndx == st_ndx);
            loo_sd = sds{loo_ndx};
            loo_tgt(iter_ndx) = loo_sd.fm10(times(loo_in_st));
            loo_tm(iter_ndx) = Ts(t);
            stt = times(loo_in_st);
            y_loo = [0.5*(loo_sd.ew(stt)+loo_sd.ed(stt)),1,loo_sd.elevation,loo_sd.rain(stt)]';

            % construct regressors, observations and variances
            if(Nst<10)
                fprintf('TSM time %g building %d regressors.\n', Ts(t), Nst-1);
            end
            X = zeros(Nst-1,4);
            Z = zeros(Nst-1,1); 
            G = zeros(Nst-1,1);
            rmap = [1:loo_in_st-1,nan,loo_in_st:Nst-1];
            for o=[1:loo_in_st-1,loo_in_st+1:Nst]
                std = sds{st_ndx(o)};
                stt = times(o);
                X(rmap(o),:) = [0.5*(std.ew(stt)+std.ed(stt)),1,std.elevation,std.rain(stt)];
                Z(rmap(o)) = std.fm10(stt);
                G(rmap(o)) = std.fm10_var(stt);
            end
            if(all(X(:,4)==0))
                X = X(:,1:3);
                y_loo = y_loo(1:3);
            end
            [beta,sigma2] = estimate_tsm_parameters(X,Z,G);
            [betac,sigma2c] = estimate_tsm_parameters_constr(X,Z,G,[X;y_loo']);
            betas(iter_ndx) = beta(1);
            betasc(iter_ndx) = betac(1);

            % find regressor vector for target station and make prediction
            loo_tsm(iter_ndx) = y_loo' * beta;
            loo_var(iter_ndx) = sigma2 + y_loo' * ((X'*diag(1./(sigma2+G))*X)\y_loo);
            loo_tsmc(iter_ndx) = y_loo' * betac;
            loo_varc(iter_ndx) = sigma2c + y_loo' * ((X'*diag(1./(sigma2c+G))*X)\y_loo);

            iter_ndx = iter_ndx + 1;
        end

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
        
        print('-dpng', [station_code,'_tsm_vs_time']);

        f2 = figure();
        set(f2,'units','centimeters');
        set(f2,'papersize',[12,12]);
        set(f2,'paperposition',[0,0,12,12]);
        
        subplot(221);
        plot(loo_tgt,loo_tsm, 'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        hold on;
        plot([0,0.5],[0,0.5],'b-');
        hold off;
        axis('square');
        title([sds{z}.stid,': Observation vs. LS TSM fit']);
        xlabel('observation [-]');
        ylabel('TSM fit [-]');

        subplot(222);
        plot(abs(loo_tgt-loo_tsm),sqrt(loo_var),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        title([sds{z}.stid,': Abs. error vs. LS TSM stdev estimate']);
        axis('square');
        xlabel('abs. error [-]');
        ylabel('stdev [-]');

        subplot(223);
        plot(loo_tgt,loo_tsmc, 'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        hold on;
        plot([0,0.5],[0,0.5],'b-');
        hold off;
        axis('square');
        title([sds{z}.stid,': Observation vs. CLS TSM fit']);
        xlabel('observation [-]');
        ylabel('TSM fit [-]');

        subplot(224);
        plot(abs(loo_tgt-loo_tsmc),sqrt(loo_varc),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        title([sds{z}.stid,': Abs. error vs. CLS TSM stdev estimate']);
        axis('square');
        xlabel('abs. error [-]');
        ylabel('stdev [-]');
        
        print('-dpng', [station_code,'_tsm_scatters']);
        
        close all;

    end
    
    % report
    Nt = iter_ndx - 1;
    fprintf('Station: %s LS/MAPE: %g  LS/MSE: %g  CLS/MAPE: %g   CLS/MSE: %g\n', ...
            station_code, norm(loo_tgt-loo_tsm,1)/Nt, norm(loo_tgt-loo_tsm,2)/Nt, ...
            norm(loo_tgt-loo_tsmc,1)/Nt, norm(loo_tgt-loo_tsmc,2)/Nt);
    
end

