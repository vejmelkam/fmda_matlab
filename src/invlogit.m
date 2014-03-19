

function T=invlogit(D,Max,N)
    E = exp(D(1:N));
    T = D;
    T(1:N) = Max*E ./ (E + 1);