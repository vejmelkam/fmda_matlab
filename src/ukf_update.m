
%
% Executes the UKF update step given the forecast and an observation.
%
%  synopsis: [ma,Pa,K,S] = ukf_update(mf,sqrtPf,H,d,R)
%
%  ARGUMENTS
%  mf     - the forecast state
%  sqrtPf - the square root of the forecast covariance
%  H      - observation operator (assumed linear)
%  d      - the data point (observation)
%  R      - the data error covariance
%
%  RETURNS
%  ma   - the assimilated state
%  Pa   - assimilated state error covariance
%  K    - the Kalman gain
%  S    - the innovation variance
%

function [ma,Pa,K,S] = ukf_update(mf,sqrtP,Pf,H,d,R)

    M = size(mf,1);
    Npts = 2*M+1;
    
    % generate new sigma points
    [sigs, w] = ukf_select_sigma_points(mf, Pf, 0);

    % generate forecast-derived observation
    y_pred_i = H * sigs;
    y_pred = sum(y_pred_i * diag(w), 2);

    % innovation covariance (H=I due to direct observation)
    sqrtS = (y_pred_i - repmat(y_pred, 1, Npts)) * diag(sqrt(w));
    S = sqrtS * sqrtS' + R; 

    % the cross covariance of state & observation errors
    Cxy = sqrtP * sqrtS';

    % Kalman gain is inv(S) * P for this case (direct observation)
    K = Cxy / S;

    % update step of Kalman filter to shift model state
    ma = mf + K*(d - y_pred);

    % state error covariance is reduced by the observation
    Pa = Pf - K*S*K';
