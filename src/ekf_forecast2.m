
%
% Computes the model forecast using the EKF propagation rule.
% Propagates the state error covariance via the Jacobian.
%
% synopsis: [mf,Pf,J,mid] = ekf_forecast2(Tk,Ed,Ew,m,r,dt,P,Qphr,S,rk,r0,Trk)
%
%    ARGUMENTS
%    Tk - time constants [hrs]
%    Ed,Ew - drying/wetting equilibria
%    m - moisture model state
%    r - the rain intensity
%    dt - time 
%    P - the covariance matrix associated with current state errors
%    Qphr - process noise incurred (per hour)
%
%    RETURNS:
%    mf     - the forecast
%    Pf     - the forecast covariance 
%    J      - the Jacobian
%    mid    - the model id
%

function [mf,Pf,J,mid] = ekf_forecast2(Tk,Ed,Ew,m,r,dt,P,Qphr,S,rk,r0,Trk,mdE)

    % forecast next point
    [mf, mid] = moisture_model_ext2(Tk,Ed,Ew,m,r,dt,S,rk,r0,Trk,mdE);
    
    % update covariance matrix using the tangent linear model
    J = moisture_tangent_model_ext(Tk,Ew,Ed,m,r,dt);
    Pf = J*P*J' + Qphr*dt/3600;
