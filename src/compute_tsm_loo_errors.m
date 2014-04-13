
%
% Evaluates the MAPE,mean error and RMSE for a leave-one-out spatial
% experiment.
%
% Synopsis: [sns,errs,aerrs,rmse] = compute_tsm_loo_errors(dirname)
%


function [sns,errs,aerrs,rmse,Nt] = compute_tsm_loo_errors(dirname,fieldname)

    S = dir([dirname,'/*.mat']);
    Nst = length(S);
    
    errs = zeros(Nst,1);
    aerrs = zeros(Nst,1);
    rmse = zeros(Nst,1);
    sns = cell(Nst,1);
    
    total = 0;
    Nt = 0;
    for i=1:Nst
        d = load([dirname,'/',S(i).name]);
        pred = getfield(d,fieldname);
        valid = isfinite(d.loo_tgt);
        N = sum(valid);
        if(N > 0)
            total = total + 1;
            sns{total} = S(i).name(1:5);
            errs(total) = sum(pred(valid) - d.loo_tgt(valid));
            aerrs(total) = norm(pred(valid) - d.loo_tgt(valid),1);
            rmse(total) = norm(pred(valid) - d.loo_tgt(valid),2);
            Nt = Nt + N;
        end
    end
    
    errs = errs(1:total);
    aerrs = aerrs(1:total);
    rmse = rmse(1:total);
    sns = sns(1:total);
    
    
    
    
    