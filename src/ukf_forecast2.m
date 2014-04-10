
%
% Computes the model forecast using the UKF propagation rule.  Re-estimates
% the forecast covariance as well.
%
% synopsis: [mf,sqrtP,Pf,mid] = ukf_forecast2(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr,S,rk,r0,Trk)
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
%
function [mf,sqrtP,Pf] = ukf_forecast2(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr,S,rk,r0,Trk)

    M = size(m,1);
    Npts = 2*M;

    m_sigma = ukf_select_sigma_points(m,P);

    f_sigma = zeros(M, Npts);
    for n=1:Npts
        f_sigma(:,n) = moisture_model_ext2(Tk,Ed,Ew,m_sigma(:,n),r,dt,dcay,S,rk,r0,Trk);
    end

    % compute the prediction mean x_mean(i|i-1)
    mf = 1/Npts*sum(f_sigma, 2);
    
    % FIXME: need to replace this by direct square root propagation
    sqrtP = (f_sigma - repmat(mf, 1, Npts));
    Pf = Qphr*dt/3600 + 1/Npts*(sqrtP * sqrtP');
