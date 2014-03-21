
%
% Computes the model forecast using the UKF propagation rule.  Re-estimates
% the forecast covariance as well.
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
%
%    RETURNS:
%    mf     - the forecast
%    sqrtP  - the square root of the forecast covariance (without process
%    noise)
%    Pf     - the forecast covariance
%    mid    - the model id

function [mf,sqrtP,Pf,mid] = ukf_forecast(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr)

    M = size(m,1);
    Npts = 2*M+1;

    [m_sigma, w] = ukf_select_sigma_points(m,P,0);

    f_sigma = zeros(M, Npts);
    for n=1:Npts-1
        f_sigma(:,n) = moisture_model_ext(Tk,Ed,Ew,m_sigma(:,n),r,dt,dcay);
    end
    [f_sigma(:,Npts),mid] = moisture_model_ext(Tk,Ed,Ew,m,r,dt,dcay);

    % compute the prediction mean x_mean(i|i-1)
    mf = sum(f_sigma * diag(w), 2);
    
    % FIXME: need to replace this by direct square root propagation
    sqrtP = (f_sigma - repmat(mf, 1, Npts)) * diag(sqrt(w));
    Pf = Qphr*dt/3600 + sqrtP * sqrtP';
