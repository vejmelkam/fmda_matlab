%
% Estimates the microscale variability variance by numerically solving
% a rational equation, see Fay&Herriot, 1979 for the model & equation.
%
% synopsis: s2 = bisect_solve(res2,gammas2,k)
%

function s2 = bisect_solve(res2,gammas2,k)

    s2l = 0;
    s2r = 0.01;
    tgt = numel(res2) - k;

    val_left = sum(res2./gammas2);
    val_right = sum(res2./(s2r + gammas2));

    % if for sigma^2 = 0 we are not above n-k, no real solutions exist
    % set estimate to zero and return
    if(val_left < tgt)
        s2 = 0;
        return;
    end
    
    % while the right value is above n-k, keep doubling it.  Eventually we
    % must get below any positive target n-k since the limit for sigma^2
    % growing toward infinity is zero.
    while(val_right > tgt)
        s2l = s2r;
        val_left = val_right;
        s2r = s2r * 2;
        val_right = sum(res2/(gammas2+s2r));
    end

    % keep halving intervals until we are close enough to the target value
    while(val_left-val_right > 1e-6)
        s2n = 0.5 * (s2l+s2r);
        val_new = sum(res2./(s2n+gammas2));
        if(val_new > tgt)
            val_left = val_new;
            s2l = s2n;
        else
            val_right = val_new;
            s2r = s2n;
        end
    end

    % return the midpoint between the limits
    s2 = 0.5 * (s2l+s2r);
end
