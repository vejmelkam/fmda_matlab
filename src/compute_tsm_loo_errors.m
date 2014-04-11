


function [sns,errs,aerrs,rmse] = compute_tsm_loo_errors(dirname)

    S = dir([dirname,'/*.mat']);
    Nst = length(S);
    
    errs = zeros(Nst,1);
    aerrs = zeros(Nst,1);
    rmse = zeros(Nst,1);
    sns = cell(Nst,1);
    
    total = 0;
    for i=1:Nst
        d = load([dirname,'/',S(i).name]);
        if(isfield(d,'loo_tsm'))
            pred = d.loo_tsm;
        else
            pred = d.loo_interp;
        end
        valid = isfinite(d.loo_tgt);
        N = sum(valid);
        if(N > 0)
            total = total + 1;
            sns{total} = S(i).name(1:5);
            errs(total) = mean(pred(valid) - d.loo_tgt(valid));
            aerrs(total) = norm(pred(valid) - d.loo_tgt(valid),1)/N;
            rmse(total) = norm(pred(valid) - d.loo_tgt(valid),2)/N;
        end
    end
    
    errs = errs(1:total);
    aerrs = aerrs(1:total);
    rmse = rmse(1:total);
    sns = sns(1:total);
    
    
    
    
    