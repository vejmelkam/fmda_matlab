
%
% Executes the UKF update step given the forecast and an observation.
%
% Implementation according to Simon,2010, pp 449-450.
%
%  synopsis: [ma,Pa,K,S] = ukf_update(mf,sqrtPf,H,d,R)
%
%  ARGUMENTS
%  mf      - the forecast state
%  sqrtPf  - the square root of the forecast covariance
%  f_sigma - sigma points with state passed through model
%  H       - observation operator (assumed linear)
%  d       - the data point (observation)
%  R       - the data error covariance
%
%  RETURNS
%  ma   - the assimilated state
%  Pa   - assimilated state error covariance
%  K    - the Kalman gain
%  S    - the innovation variance
%

function [ma,Pa,K,S] = ukf_update(mf,sqrtP,sigma_f,H,d,R,kappa)

    nx = size(sqrtP,1);
    n = size(sigma_f,1);
    nv = size(H,1);
    Npts = 2*n+1;
    
    % compute the weights for given kappa
    w = ones(Npts,1) * 1/(2*(n+kappa));
    w(Npts) = kappa / (n+kappa);
    
    % the sigma_f points are missing the noise component, add it now
    sigma_f(end-nv+1:end,n) = sqrt(n+kappa)*sqrt(R);
    sigma_f(end-nv+1:end,2*n) = -sqrt(n+kappa)*sqrt(R);

    % generate forecast-derived observation
    y_pred_i = H * sigma_f(1:nx,:) + sigma_f(end-nv+1:end,:);
    y_pred = y_pred_i*w;

    % innovation covariance
    sqrtS = (y_pred_i - y_pred)*diag(w.^0.5);
    S = (sqrtS * sqrtS'); 

    % the cross covariance of state & observation errors
    Cxy = sqrtP * sqrtS';

    % Kalman gain is inv(S) * P for this case (direct observation)
    K = Cxy / S;

    % update step of Kalman filter to shift model state
    ma = mf(1:nx) + K*(d - y_pred);

    % state error covariance is reduced by the observation
    Pa = (sqrtP*sqrtP') - K*S*K';
