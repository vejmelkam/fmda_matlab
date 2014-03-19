

function T = logit(D,Max,N)
    T = D;
    T(1:N) = log(D(1:N) ./ (Max - D(1:N)));