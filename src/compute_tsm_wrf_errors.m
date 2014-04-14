
%
% Evaluates the MAPE,mean error and RMSE for the WRF forecast test.
%
% Synopsis: [sns,errs,aerrs,rmse] = compute_tsm_wrf_errors(dirname,stop_da_at)
%

function [sns,errs,aerrs,rmse,Nt] = compute_tsm_wrf_errors(fname,rng)

    d = load(fname);
    Nst = length(d.sds);
    
    errs = zeros(Nst,1);
    aerrs = zeros(Nst,1);
    rmse = zeros(Nst,1);
    sns = cell(Nst,1);
    Nt = zeros(Nst,1);
    
    for i=1:Nst
        fm10 = d.fm10_tgt(:,i);
        pred = d.fm10_model(:,i,2);
        valid = false(size(fm10));
        valid(rng) = true;
        valid(~isfinite(fm10)) = false;
        N = sum(valid);
        if(N > 0)
            Nt(i) = N;
            sns{i} = d.sds{i}.stid;
            errs(i) = sum(pred(valid) - fm10(valid));
            aerrs(i) = norm(pred(valid) - fm10(valid),1);
            rmse(i) = norm(pred(valid) - fm10(valid),2);
        end
    end
        
    
    
    
    