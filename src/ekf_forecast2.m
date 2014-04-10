
%
% Computes the model forecast using the EKF propagation rule.
% Propagates the state error covariance via the Jacobian.
%
% synopsis: [mf,Pf,J,mid] = ekf_forecast2(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr,S,rk,r0,Trk)
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
%    Pf     - the forecast covariance 
%    J      - the Jacobian
%    mid    - the model id
%

function [mf,Pf,J,mid] = ekf_forecast2(Tk,Ed,Ew,m,r,dt,dcay,P,Qphr,S,rk,r0,Trk)

    % forecast next point
    [mf, mid] = moisture_model_ext2(Tk,Ed,Ew,m,r,dt,dcay,S,rk,r0,Trk);
    
    % update covariance matrix using the tangent linear model
    J = moisture_tangent_model_ext(Tk,Ew,Ed,m,r,dt,dcay);
    Pf = J*P*J' + Qphr*dt/3600;
