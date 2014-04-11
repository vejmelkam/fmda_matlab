
%
% Computes the model forecast using the UKF propagation rule.  Re-estimates
% the forecast covariance as well.
%
% Implementation according to Simon,2010, pp 449-450.
%
% synopsis: [mf,sqrtP,Pf,mid] = ukf_forecast(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr)
%
%    ARGUMENTS
%    Tk - time constants [hrs]
%    Ed,Ew - drying/wetting equilibria
%    m - moisture model state
%    r - the rain intensity
%    dt - time 
%    dcay - the decay time constant [hrs] of the assimilated coefficients
%    P - the covariance matrix associated with current state errors
%    Qphr - process noise incurred (per hour)
%    R - the observation noise covariance
%    kappa - variable allowing incorporation of prior information into the
%            Unscented transform, set to 0 if in doubt
%
%    RETURNS:
%    mf     - the forecast
%    sqrtP  - the square root of the forecast covariance (without process
%    noise)
%    Pf     - the forecast covariance

function [mf,sqrtP,f_sigma] = ukf_forecast(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr,R,kappa)

    nx = size(m,1);
    nv = size(R,1);
    n = 2*nx+ nv;
    Npts = 2*n+1;

    % construct sigma points from last assimilated state
    f_sigma = ukf_select_sigma_points(m,P,Qphr*dt/3600,R,kappa);

    % build weights for the sigma points
    w = ones(Npts,1) * 1/(2*(n+kappa));
    w(Npts) = kappa / (n+kappa);

    % pass sigma points through the model function (only the state!)
    % noise terms remain stored as constructed in f_sigma above
    for i=1:Npts
        f_sigma(1:nx,i) = moisture_model_ext(Tk,Ed,Ew,f_sigma(1:nx,i),r,dt,dcay) + f_sigma(nx+1:2*nx,i);
    end

    % compute the prediction mean x_mean(i|i-1)
    mf = f_sigma * w;
    
    % compute sqryP*sqrtP' = Pf (forecast covariance)
    sqrtP = (f_sigma(1:nx,:) - repmat(mf(1:nx), 1, Npts))*diag(w.^0.5);
