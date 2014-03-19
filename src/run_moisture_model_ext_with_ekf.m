
%
%  This script runs the extended moisture model for one grid point with the
%  Extended Kalman filter.
%

% Load station data & compute equilibria
station = 'BAWC2';
yr = '2012';
sd = load_station_data(station,yr,[2400,4800]);
tdays = sd.tdays - sd.tdays(1);
Ed = sd.ed;
Ew = sd.ew;

% time in hours (integration step for moisture is 1h)
N = length(tdays);
do_assim = false(N,1);
do_assim(1:12:N) = true;

% number of fuels and extended model size
k = 3;
M = k+2;

% parameters of the simulation
m_ext = zeros(M,1);
m_ext(1:k) = 0.5 * (Ed(1)+Ew(1));

P = eye(M) * 0.01;   % error covariance of the initial guess

% Kalman filter Q (model error covariance)
Qphr = zeros(M);
Qphr(1:k, 1:k) = diag([0.0005,0.00005,0.00001]);
Qphr(k+1,k+1) = 0.0001;
Qphr(k+2,k+2) = 0.0001;

% the observation operator is a n_k x M matrix with I_(n_k) on the left
H = zeros(1,M);
H(1,2) = 1;

% storage space for results (with filtering)
m_f = zeros(N, M);
m_f(1, :) = m_ext';

m_n = zeros(N, M); % (without filtering)
m_n(1, :) = m_ext';

m_c = zeros(N, M); % (without filtering, with corrected dE, dS)
m_c(1, :) = m_ext';
m_c(1, k+1) = -0.05;
m_c(1, k+2) = -0.6;

% indicator of moisture model that is switched on at each time point
model_ids = zeros(N, k);

% storage for matrix traces
trP = zeros(N, 1);
sK = zeros(N, M);
trS = zeros(N, 1);
trJ = zeros(N, 1);
sP = zeros(N, M, M);
sJ = zeros(N, M, M);

% predict & update 0loop
for i=2:N
    
    % compute the integration time step in seconds
    dt = (tdays(i) - tdays(i-1)) * 86400;
    Ed2 = (Ed(i)+Ed(i-1))/2;
    Ew2 = (Ew(i)+Ew(i-1))/2;
    
    % compute & store results for system without Kalman filtering
    m_n(i, :) = moisture_model_ext([1,10,100]',Ed2,Ew2,m_n(i-1,:)',sd.rain(i),dt,1e10);
    
    % compute & store results for system without Kalman filtering but with
    % initial dE & dS correction
    m_c(i, :) = moisture_model_ext([1,10,100]',Ed2,Ew2,m_c(i-1,:)',sd.rain(i),dt,1e10);

    % KALMAN PREDICT STEP
    
    % estimate new moisture mean based on last best guess (m)
    [m_pred, model_ids(i,:)] = moisture_model_ext([1,10,100]',Ed2,Ew2,m_f(i-1,:)',sd.rain(i),dt,1e10);
    
    % update covariance matrix using the tangent linear model
    Jm = moisture_tangent_model_ext([1,10,100]',Ew2,Ed2,m_f(i-1,:)',sd.rain(i),dt,1e10);
    sJ(i, :, :) = Jm;
    P = Jm*P*Jm' + Qphr*dt/3600;
    sP(i, :, :) = P;
    %trP(i) = trace(P);
    trP(i) = prod(eig(P(1:k,1:k)));
    trJ(i) = prod(eig(Jm));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if(do_assim(i))
        
        % acquire current measurement & move to next one
        m_measured = sd.fm10(i);
        R = sd.fm10_var(i);

        % innovation covariance
        S = H*P*H' + R;
        trS(i) = prod(eig(S));
        
        % Kalman gain is inv(S) * P for this case (direct observation)
        K = P * H' / S;
        sK(i,:) = K;
        %trK(i) = prod(eig(K(1:k,1:k)));
        
        % update step of Kalman filter to shift model state
        m_f(i,:) = m_pred + K*(m_measured - H*m_pred);
        
        % state error covariance is reduced by the observation
        P_red = K*S*K';
        P = P - P_red;
    
    else
        
        % if no observation is available, store the predicted value
        m_f(i,:) = m_pred;
        
    end
        
end

set(0,'DefaultAxesLooseInset',[0 0 0 0])

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(311);
plot(tdays, m_f(:,2), 'r-', 'linewidth', 1.5);
hold on;
plot(tdays, m_n(:,2), 'g-', 'linewidth', 1.5);
plot(tdays, m_c(:,2), 'b-', 'linewidth', 1.5);
%plot(tdays, rain, 'k--', 'linewidth', 1.5);
plot(tdays(do_assim), sd.fm10(do_assim), 'ko', 'markersize', 4, 'markerfacecolor', 'm');
%plot(repmat(tdays, 1, 2), [m_f(:,2) - sqrt(sP(:, 2, 2)), m_f(:,2) + sqrt(sP(:, 2, 2))], 'rx');
h = legend('sys + EKF', 'sys', 'sys + corr', 'obs', 'orientation', 'horizontal');
set(h, 'fontsize', 10);
title('Plot of the evolution of the moisture model [EKF]', 'fontsize', 12);
ylim([0, min(1.2,1.1*max([m_f(:,2);m_n(:,2);m_c(:,2)]))]);

% select time indices corresponding to observation times
subplot(312);
plot(tdays, log10(trP), 'b-', 'linewidth', 1.5);
hold on;
plot(tdays, log10(trJ), 'k-', 'linewidth', 1.5);
plot(tdays(do_assim), log10(trS(do_assim)), 'ko', 'markerfacecolor', 'green', 'markersize', 4);
plot(tdays(do_assim), log10(abs(sK(do_assim,2))), 'ko', 'markerfacecolor', 'red', 'markersize', 4);
plot(tdays(do_assim), log10(abs(sK(do_assim,1))), 'ko', 'markerfacecolor', 'blue', 'markersize', 4);
plot(tdays(do_assim), log10(abs(sK(do_assim,3))), 'ko', 'markerfacecolor', 'cyan', 'markersize', 4);
hold off;
a = axis();
axis([a(1) a(2) a(3) max(a(4), 1.0)]);
h = legend('State', 'Jacobian', 'Innovation', 'K(2)', 'K(1)', 'K(3)', 'orientation', 'horizontal');
set(h, 'fontsize', 10);
title('Kalman filter: log(generalized variance) of covariances & K(2) element  time [EKF]', 'fontsize', 12);
ylim([-12, 5]);

subplot(313);
plot(tdays, sP(:, 2, 2), 'r-', 'linewidth', 1.5);
hold on
plot(tdays, sP(:, 2, 4), 'g-', 'linewidth', 1.5);
plot(tdays, sP(:, 2, 5), 'b-', 'linewidth', 1.5);
plot(tdays, sP(:, 1, 1), 'k-', 'linewidth', 1.5);
plot(tdays, sP(:, 3, 3), 'c-', 'linewidth', 1.5);
hold off
h = legend('var(m(2))', 'cov(m,dE)', 'cov(m,dS)', 'var(m(1))', 'var(m(3))', 'orientation', 'horizontal');
set(h, 'fontsize', 12);
title('Covariance between moisture and system parameters [EKF]', 'fontsize', 12);

%print(gcf, '-depsc', 'ekf_assim_ts.eps');

figure('units','normalized','outerposition',[0 0 1 1])
subplot(311);
plot(repmat(tdays, 1, 2), m_f(:, [4,5]), 'linewidth', 1.5);
h = legend('dE', 'dS');
set(h, 'fontsize', 12);
title('Equilibrium changes [EKF]', 'fontsize', 12);

subplot(312);
plot(repmat(tdays, 1, 2), [Ed + m_f(:,4),Ew + m_f(:,4)], 'linewidth', 1.5);
h = legend('Ed+dE', 'Ew+dE');
set(h, 'fontsize', 12);
title('Adjusted equilibria [EKF]', 'fontsize', 12);

subplot(313);
plot(tdays, model_ids(:,2), 'or', 'markerfacecolor', 'red');
title('Active submodel of the moisture model for 10-hr fuel [EKF]', 'fontsize', 12);
ylim([-1, 5]);

%print(gcf, '-depsc', 'ekf_assim_params.eps');

fprintf('MSE for EKF: %g  MAPE for EKF: %g\n', norm(sd.fm10 - m_f(:,2),2)/N,norm(sd.fm10 - m_f(:,2),1)/N);
