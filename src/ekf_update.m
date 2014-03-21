%
% Executes the UKF update step given the forecast and an observation.
%
%  synopsis: [ma,Pa,K,S] = ekf_update(mf,Pf,H,d,R)
%
%  ARGUMENTS
%  mf     - the forecast state
%  H      - observation operator (assumed linear)
%  d      - the data point (observation)
%  R      - the data error covariance
%
%  RETURNS
%  ma   - the assimilated state
%  P    - assimilated state error covariance
%  K    - the Kalman gain
%  S    - the innovation variance
%

function [ma,Pa,K,S] = ekf_update(mf,Pf,H,d,R)

    % innovation covariance
    S = H*Pf*H' + R;
        
    % Kalman gain is inv(S) * P for this case (direct observation)
    K = Pf * H' / S;
        
    % update step of Kalman filter to shift model state
    ma = mf + K*(d - H*mf);
        
    % state error covariance is reduced by the observation
    Pa = Pf - K*S*K';