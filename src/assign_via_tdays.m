
%
% Given an observation type x1 that is available at times t1,
% assign this to a variable y available at times t2.  When no observation
% is available, the value of y is nan, otherwise the appropriate value from
% x1 is compied.
%
% synopsis: y = assign_via_tdays(t1,x1,t2)
%
%

function y = assign_via_tdays(t1,x1,t2)
    N1 = length(t1);
    N2 = length(t2);
    y = nan * zeros(N2,1);
    
    i1 = 1;
    i2 = 1;
    while((i1 <= N1) && (i2 <= N2))
        if(t1(i1) == t2(i2))
            y(i2) = x1(i1);
            i1 = i1 + 1;
            i2 = i2 + 1;
        elseif(t1(i1) < t2(i2))
            i1 = i1 + 1;
        else
            i2 = i2 + 1;
        end
    end
    