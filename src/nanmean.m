

function m = nanmean(X)

    m = zeros(1,size(X,2));
    for i=1:size(X,2)
        mask = isfinite(X(:,i));
        m(i) = sum(X(mask,i)) / sum(mask);
    end