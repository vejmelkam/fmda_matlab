


A = [0.3,0.5;0.3,-0.2];
x0 = [1;1];
qd = [0.4;0.2];
Q = diag(qd);
N = 50;
R = 0.3;
kappa = 1;
Npts = 11;

H = [1,0];

xt = zeros(N,2);
xu = zeros(N,2);
xe = zeros(N,2);

xt(1,:) = x0;
xu(1,:) = x0;
xu2(1,:) = x0;
xe(1,:) = x0;

P0 = diag([0.01,0.01]);

Pe = zeros(N,2,2);
Pu = zeros(N,2,2);

Pe(1,:,:) = P0;
Pu(1,:,:) = P0;

Kes = zeros(N,2);
Kus = zeros(N,2);

Ses = zeros(N,1);
Sus = zeros(N,1);

Au = zeros(2,5);
Au(1:2,1:2) = A;
Au(1:2,3:4) = eye(2);

for i=2:N
    
    % move true model one timestep
    eta = randn(2,1) .* sqrt(qd);
    xt(i,:) = A*xt(i-1,:)' + eta;
    
    % construct 'observation'
    d = H*xt(i,:)' + randn()*sqrt(R);
    
    % forecast via UKF
    [mu,sqrtPu,sigma_f] = ukf_forecast_general(xu(i-1,:)',@(x,w) A*x + w,squeeze(Pu(i-1,:,:)),Q,1,kappa);

%     % propagate via UKF
%     Pui(1:2,1:2) = Pu(i-1,:,:);
%     m_sigma  = ukf_select_sigma_points(xu(i-1,:)',squeeze(Pu(i-1,:,:)),Q,R,kappa);
% 
%     w = ones(Npts,1) * 1/(2*(5+kappa));
%     w(Npts) = kappa / (5+kappa);
%     
%     f_u = Au*m_sigma;
% 
%     % compute the prediction mean x_mean(i|i-1)
%     xu(i,:) = sum(f_u * diag(w), 2);
%     
%     sqrtP = (f_u - repmat(xu(i,:)', 1, Npts))*diag(w.^0.5);
%     Pf = (sqrtP * sqrtP');
    xu(i,:) = mu;
    Pu(i,:,:) = sqrtPu*sqrtPu';
    
    % propagate via EKF
    Pei = A*squeeze(Pe(i-1,:,:))*A' + Q;
    Pe(i,:,:) = Pei;
    xe(i,:) = A*xe(i-1,:)';
    
%     % generate forecast-derived observation
%     y_pred_i = H * f_u + m_sigma(5,:);
%     y_pred = y_pred_i * w;
% 
%     % innovation covariance
%     sqrtS = (y_pred_i - y_pred)*diag(w.^0.5);
%     Su = (sqrtS * sqrtS');
%     Sus(i) = Su;
% 
%     % the cross covariance of state & observation errors
%     Cxy = sqrtP * sqrtS';
% 
%     % Kalman gain is inv(S) * P for this case (direct observation)
%     Ku = Cxy / Su;
%     Kus(i,:) = Ku;
% 
%     % update step of Kalman filter to shift model state
%     xu(i,:) = xu(i,:)' + Ku*(d - y_pred);
% 
%     % state error covariance is reduced by the observation
%     Pu(i,:,:) = squeeze(Pu(i,:,:)) - Ku*Su*Ku';

    [ma,Pa,Ku,Su] = ukf_update(mu,sqrtPu,sigma_f,H,d,R,kappa);
    xu(i,:) = ma;
    Pu(i,:,:) = Pa;
    Kus(i,:) = Ku;
    Sus(i) = Su;
    
    % assimilate via EKF
    Pe_i = squeeze(Pe(i,:,:));
    Se = (H*Pe_i*H' + R);
    Ses(i) = Se;
    Ke = Pe_i*H'/Se;
    xe(i,:) = xe(i,:)' + Ke*(d - H*xe(i,:)');
    Pe(i,:,:) = Pe_i - Ke*Se*Ke';
    Kes(i,:) = Ke;
    
end


t = 1:N;
figure;
subplot(1,2,1);
plot(t,xt(:,1),'r-',t,xe(:,1),'g-',t,xu(:,1),'b-');
legend('ground truth','EKF','UKF');
title('First component');
subplot(1,2,2);
plot(t,xt(:,2),'r-',t,xe(:,2),'g-',t,xu(:,2),'b-');
legend('ground truth','EKF','UKF');
title('Second component');


figure;
subplot(1,3,1);
plot(t,squeeze(Pe(:,1,1)),'g-',t,squeeze(Pu(:,1,1)),'b-');
legend('EKF','UKF');
title('P(1,1) vs. time');
subplot(1,3,2);
plot(t,Kes(:,1),'g-',t,Kus(:,1),'b-');
legend('EKF','UKF');
title('K(1,1) vs. time');
subplot(1,3,3);
plot(t,Ses(:,1),'g-',t,Sus(:,1),'b-');
legend('EKF','UKF');
title('Innov. covariance vs. time');
