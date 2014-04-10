

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
    loo_ndx    = z;
    loo_tgt    = zeros(numel(Ts),1);
    loo_interp = zeros(numel(Ts),1);
    loo_tm     = zeros(numel(Ts),1);
    loo_var    = zeros(numel(Ts),1);
    loo_sd = sds{z};
    station_code = sds{z}.stid;
    fprintf('Leaving out station %d: %s\n', z, station_code);
    
    iter_ndx = 1;
    for t=1:length(Ts)

        [ndxs,fm10o] = find_valid_obs(Ts(t),sds);
        ndxs(isnan(fm10o)) = nan;
        st_ndx = find(isfinite(ndxs));
        Nst = length(st_ndx);
        sts_per_time(t) = Nst;
        times = ndxs(st_ndx);

        if((length(st_ndx) > 1) && any(st_ndx == loo_ndx))
            % find our left-out station in the array
            loo_in_st = find(loo_ndx == st_ndx);
            loo_sd = sds{loo_ndx};
            loo_tgt(iter_ndx) = loo_sd.fm10(times(loo_in_st));
            loo_tm(iter_ndx) = Ts(t);
            stt = times(loo_in_st);

            % construct regressors, observations and variances
            D = zeros(Nst-1,1);
            Z = zeros(Nst-1,1); 
            G = zeros(Nst-1,1); 
            rmap = [1:loo_in_st-1,nan,loo_in_st:Nst-1];
            for o=[1:loo_in_st-1,loo_in_st+1:Nst]
                std = sds{st_ndx(o)};
                stt = times(o);
                D(rmap(o),:) = great_circle_distance(loo_sd.lat,loo_sd.lon,std.lat,std.lon);
                Z(rmap(o)) = std.fm10(stt);
                G(rmap(o)) = std.fm10_var(stt);
            end
            
            % compute weights, interpolation and variance
            Id2 = D.^-2;
            W = Id2 / sum(Id2);
            loo_interp(iter_ndx) = W' * Z;
            loo_var(iter_ndx) = sum(W.^2 .* G);

            iter_ndx = iter_ndx + 1;
        end

    end
    
    loo_tgt    = loo_tgt(1:iter_ndx-1);
    loo_interp = loo_interp(1:iter_ndx-1);
    loo_var    = loo_var(1:iter_ndx-1);
    loo_tm     = loo_tm(1:iter_ndx-1);
    
    if(iter_ndx > 5)
    
        f1 = figure();
        set(f1,'units','centimeters');
        set(f1,'papersize',[18,6])        
        set(f1,'paperposition',[0,0,18,6]);
        
        subplot(211);
        plot([loo_tm,loo_tm],[loo_tgt,loo_interp],'linewidth',1.2);
        legend('tgt', 'interp');
        title([sds{z}.stid,': Comparison between TSM fit (LS and CLS) and observation']);
        xlabel('Time [days from 1.1]');
        ylabel('fm10 [-]');

        subplot(212);
        plot([loo_tm,loo_tm],[abs(loo_tgt-loo_interp),sqrt(loo_var)]);
        title([sds{z}.stid,': INTERP2-obs error vs. INTERP variance estimate']);
        legend('abs err', 'stdev');
        xlabel('Time [days from 1.1]');
        ylabel('fm10 [-]');
        
        print('-dpng', [station_code,'_interp2_vs_time']);

        f2 = figure();
        set(f2,'units','centimeters');
        set(f2,'papersize',[12,12]);
        set(f2,'paperposition',[0,0,12,12]);
        
        subplot(121);
        plot(loo_tgt,loo_interp, 'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        hold on;
        plot([0,0.5],[0,0.5],'b-');
        hold off;
        axis('square');
        title([sds{z}.stid,': Observation vs. INTERP2 fit']);
        xlabel('observation [-]');
        ylabel('TSM fit [-]');

        subplot(122);
        plot(abs(loo_tgt-loo_interp),sqrt(loo_var),'o', 'MarkerFaceColor',[0,0,0],'MarkerSize',5)
        title([sds{z}.stid,': Abs. error vs. INTERP2 stdev estimate']);
        axis('square');
        xlabel('abs. error [-]');
        ylabel('stdev [-]');
        
        print('-dpng', [station_code,'_interp_scatters']);
        
        % report
        valids = find(isfinite(loo_tgt));
        Nt = length(valids);
        fprintf('Station: %s INTERP2/MAPE: %g  INTERP2/MSE: %g\n', ...
                station_code, norm(loo_tgt(valids)-loo_interp(valids),1)/Nt, ...
                norm(loo_tgt(valids)-loo_interp(valids),2)/Nt);

        % dump the results to file
        sdsz = sds{z};
        save([station_code,'.mat'],'loo_tm','loo_interp','loo_tgt','sdsz','loo_var');
        
    else
        fprintf('Not enough valid points for station %s\n', station_code);
    end
    
    
end

