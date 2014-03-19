

function ndxs = find_valid_obs(t,sds)

    N = length(sds);
    ndxs = zeros(N,1);
    
    for i=1:N
        n = find(abs(sds{i}.tdays-t) < 1/48);
        if(isempty(n))
            ndxs(i) = nan;
        else
            ndxs(i) = n(1);
        end
    end
    