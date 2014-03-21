

function [ndxs,fm10o,fm10v] = find_valid_obs(t,sds)

    N = length(sds);
    ndxs = zeros(N,1);
    fm10o = zeros(N,1);
    fm10v = zeros(N,1);
    
    for i=1:N
        [dist,n] = min(abs(sds{i}.tdays-t));
        if(dist < 1/47)
            ndxs(i) = n;
            fm10o(i) = sds{i}.fm10(n);
            fm10v(i) = sds{i}.fm10_var(n);
        else
            ndxs(i) = nan;
            fm10o(i) = nan;
            fm10v(i) = nan;
        end
    end
    