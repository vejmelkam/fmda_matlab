
%
% Computes the model forecast using the UKF propagation rule.  Re-estimates
% the forecast covariance as well.
%
% Implementation according to Simon,2010, pp 449-450.
%
% synopsis: [mf,sqrtP,f_u] = ukf_forecast(xa,f(x,w),P,Q,R,kappa)
%
%    ARGUMENTS
%    xa       - previous analysis
%    f(x,w)   - function handle that computes forecast from state & proc noise
%    P        - state error covariance
%    Q        - process noise covarance (must be diagonal)
%    R        - observation noise covariance (must be diagonal)
%    kappa - variable allowing incorporation of prior information into the
%            Unscented transform, set to 0 if in doubt
%
%    RETURNS:
%    mf      - the forecast
%    sqrtP   - the square root of the forecast covariance (without process
%    noise)
%    f_u     - sigma points passed through
%

function [xf,sqrtP,sigma_f] = ukf_forecast_general(xa,f,P,Q,nv,kappa)

    nx = size(P,1);
    n = 2*nx+ nv;
    Npts = 2*n+1;

    % construct sigma points from last assimilated state
    sigma_f = ukf_select_sigma_points(xa,P,Q,nv,kappa);

    % build weights for the sigma points
    w = ones(Npts,1) * 1/(2*(n+kappa));
    w(Npts) = kappa / (n+kappa);

    % pass sigma points through the model function (only the state!)
    % noise terms remain stored as constructed in f_sigma above
    for i=1:Npts
        sigma_f(1:nx,i) = f(sigma_f(1:nx,i),sigma_f(nx+1:2*nx,i));
    end

    % compute the prediction mean x_mean(i|i-1)
    xf = sigma_f(1:nx,:) * w;
    
    % compute sqrtP*sqrtP' = Pf (forecast covariance)
    sqrtP = (sigma_f(1:nx,:) - repmat(xf, 1, Npts))*diag(w.^0.5);
