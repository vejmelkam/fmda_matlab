
%
% Evaluates the MAPE,mean error and RMSE for the WRF forecast test.
%
% Synopsis: [sns,errs,aerrs,rmse] = compute_tsm_wrf_errors(fname,rng)
%

function [sns,errs,aerrs,rmse,Nt] = compute_tsm_wrf_errors(fname,rng,bfix)

    if(nargin < 3)
        bfix = 0;
    end

    d = load(fname);
    Nst = length(d.sds);
    Nr = length(rng);
    
    errs = zeros(Nr,Nst);
    aerrs = zeros(Nr,Nst);
    rmse = zeros(Nr,Nst);
    sns = cell(Nst);
    Nt = zeros(Nr,Nst);
    
    for i=1:Nst
        sns{i} = d.sds{i}.stid;
        fm10 = d.fm10_tgt(:,i);
        pred = d.fm10_model(:,i,2) + bfix;
        for t=1:Nr
            r = rng(t);
            if(isfinite(fm10(r)))
                Nt(t,i) = Nt(t,i)+1;
                errs(t,i) = pred(r) - fm10(r);
                aerrs(t,i) = abs(pred(r) - fm10(r));
                rmse(t,i) = (pred(r) - fm10(r))^2;
            end
        end
    end
    
