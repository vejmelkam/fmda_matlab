

%
%  This script runs the extended moisture model for one grid point with the
%  Unscented Kalman Filter.
%

% Load station data & compute equilibria
station = 'ESPC2';
yr = '2012';
[tdays,t2,relh,fm10,fm10_var,rain,accp] = load_station_data(station,yr,[1,1000]);
[Ed,Ew] = equilibrium_moisture2(relh,t2);
tdays = tdays - tdays(1);

% time in hours (integration step for moisture is 1h)
N = length(tdays);
do_assim = false(N,1);
do_assim(1:48:N) = true;

% number of fuels and extended model size
k = 3;
M = k+2;

% parameters of the simulation
m_ext = zeros(M,1);
m_ext(1:k) = 0.5 * (Ed(1)+Ew(1));

% initialize covariance
P = eye(M) * 0.01;   % error covariance of the initial guess

% Kalman filter Q (model error covariance)
Qphr = zeros(M);
Qphr(1:k, 1:k) = diag([0.0005,0.0001,0.00005]) * 2;
Qphr(k+1,k+1) = 0.0005;
Qphr(k+2,k+2) = 0.0005;

% the observation operator is a n_k x Ndim matrix with I_(n_k) on the left
H = zeros(1,k+2);
H(2) = 1;

% storage space for results (with filtering)
m_f = zeros(N, M);
m_f(1, :) = m_ext';
m_f(1, k+1) = -0.05;
m_f(1, k+2) = -0.6;

m_n = zeros(N, M); % (without filtering)
m_n(1, :) = m_ext;

m_c = zeros(N, M); % (without filtering, with initial correction)
m_c(1, :) = m_ext;
m_c(1, 4) = -0.05;
m_c(1, 5) = -0.6;

% indicator of moisture model that is switched on at each time point
model_ids = zeros(N, k);

% storage for matrix traces
trP = zeros(N, 1);
trS = zeros(N, 1);
sK = zeros(N, M);
sP = zeros(N, M, M);

% W0 is a UKF parameter affecting the sigma point distribution
Npts = M * 2 + 1;
W0 = 0;

% Storage of the sigma points
sigma_pts = zeros(N, M, Npts);

% predict & update loop
for i = 2:N
    
    % compute the integration time step
    dt = (tdays(i) - tdays(i-1)) * 86400;
    
    % draw the sigma points
    [m_sigma_lgt, w] = ukf_select_sigma_points(logit(m_f(i-1,:)',2.5,k), P, W0);
    sigma_pts(i, :, :) = m_sigma_lgt;
    
    % compute & store results for system without KF
    m_n(i, :) = moisture_model_ext([1,10,100]',Ed(i),Ew(i),m_n(i-1,:)',rain(i),dt,1e10);
    
    % compute & store results for corrected system without KF
    m_c(i, :) = moisture_model_ext([1,10,100]',Ed(i),Ew(i),m_c(i-1,:)',rain(i),dt,1e10);

    % UKF prediction step - run the sigma set through the nonlinear
    % function
    
    % estimate new moisture mean based on last best guess (m)
    m_sigma_lgt1 = zeros(M, Npts);
    for n=1:Npts-1
        m_sigma_lgt1(:,n) = logit(moisture_model_ext([1,10,100]',Ed(i),Ew(i),invlogit(m_sigma_lgt(:,n),2.5,k),rain(i),dt,1e10),2.5,k);
    end
    [m_sigma_lgt1(:,Npts), model_ids(i,:)] = moisture_model_ext([1,10,100]',Ed(i),Ew(i),m_f(i-1,:)',rain(i),dt,1e10);
    m_sigma_lgt1(:,Npts) = logit(m_sigma_lgt1(:,Npts),2.5,k);
    
    % compute the prediction mean x_mean(i|i-1)
    m_pred_lgt = sum(m_sigma_lgt1 * diag(w), 2);
    
    % estimate covariance matrix using the sigma point set
    sqrtP = (m_sigma_lgt1 - repmat(m_pred_lgt, 1, Npts)) * diag(sqrt(w));
    P = Qphr*dt/3600 + sqrtP * sqrtP';
    sP(i, :, :) = P;
    trP(i) = prod(eig(P));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if(do_assim(i))
        
        % acquire current measurement
        m_obs_lgt = logit(fm10(i),2.5,1);
        R = fm10_var(i)*500;

        % generate forecast-derived observation (in logit space)
        y_pred_lgt_i = H * m_sigma_lgt1;
        y_pred_lgt = sum(y_pred_lgt_i * diag(w), 2);
        
        % innovation covariance (in logit space)
        sqrtS = (y_pred_lgt_i - repmat(y_pred_lgt, 1, Npts)) * diag(sqrt(w));
        S = sqrtS * sqrtS' + R; 
        trS(i) = prod(eig(S));
        
        % the cross covariance of state & observation errors (logit space)
        Cxy = (m_sigma_lgt1 - repmat(m_pred_lgt, 1, Npts)) * diag(w) * (y_pred_lgt_i - repmat(y_pred_lgt, 1, Npts))';
        
        % Kalman gain is inv(S) * P for this case (direct observation)
        K = Cxy / S;
        
        %trK(i) = trace(K);
        sK(i,:) = K;
        
        % update step of Kalman filter to shift model state
        m_f(i,:) = invlogit(m_pred_lgt + K*(m_obs_lgt - y_pred_lgt),2.5,k);
        
        % state error covariance is reduced by the observation (everything
        % in logit space)
        P = P - K*S*K';
    
        % replace the stored covariance by the updated covariance after
        % processing the measurement
        trP(i) = prod(eig(P));
        
    else
        
        % if no observation is available, store the predicted value
        m_f(i,:) = invlogit(m_pred_lgt,2.5,k);
        
    end
        
end

set(0,'DefaultAxesLooseInset',[0,0,0,0])

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(311);
plot(tdays, m_f(:,2), 'r-', 'linewidth', 1.5);
hold on;
plot(tdays, m_n(:,2), 'g-', 'linewidth', 1.5);
plot(tdays, m_c(:,2), 'b-', 'linewidth', 1.5);
%plot(tdays, rain, 'k--', 'linewidth', 1.5);
plot(tdays(do_assim), fm10(do_assim), 'ko', 'markersize', 4, 'markerfacecolor', 'm');
%plot(repmat(tdays, 1, 2), [m_f(:,2) - sqrt(sP(:, 2, 2)), m_f(:,2) + sqrt(sP(:, 2, 2))], 'rx');
h = legend('sys + EKF', 'sys', 'sys + corr', 'obs', 'orientation', 'horizontal');
set(h, 'fontsize', 10);
title('Plot of the evolution of the moisture model [UKF]', 'fontsize', 12);
ylim([0, min(1.2,1.1*max([m_f(:,2);m_n(:,2);m_c(:,2)]))]);

% select time indices corresponding to observation times
subplot(312);
plot(tdays, log10(trP), 'b-', 'linewidth', 2);
hold on;
plot(tdays(do_assim), log10(trS(do_assim)), 'ko', 'markerfacecolor', 'green', 'markersize', 4);
plot(tdays(do_assim), log10(sK(do_assim,2)), 'ko', 'markerfacecolor', 'red', 'markersize', 4);
plot(tdays(do_assim), log10(sK(do_assim,1)), 'ko', 'markerfacecolor', 'blue', 'markersize', 4);
plot(tdays(do_assim), log10(sK(do_assim,3)), 'ko', 'markerfacecolor', 'cyan', 'markersize', 4);
hold off;
h = legend('state', 'innov.', 'K(2)', 'K(1)', 'K(3)', 'orientation', 'horizontal');
set(h, 'fontsize', 10);
title('Kalman filter: log(generalized variance) of covar/Kalman matrices vs. time [UKF]', 'fontsize', 12);
ylim([-12, 5]);

subplot(313);
plot(tdays, sP(:, 2, 2), 'r-', 'linewidth', 2);
hold on
plot(tdays, sP(:, 2, 4), 'g-', 'linewidth', 2);
plot(tdays, sP(:, 2, 5), 'b-', 'linewidth', 2);
plot(tdays, sP(:, 1, 1), 'k-', 'linewidth', 2);
plot(tdays, sP(:, 3, 3), 'c--', 'linewidth', 2);
hold off
h = legend('var(m(2))', 'cov(m(2),dE)', 'cov(m(2),dS)', 'var(m(1))', 'var(m(3))', 'orientation', 'horizontal');
set(h, 'fontsize', 12);
title('Covariance between moisture and system parameters [UKF]', 'fontsize', 12);
ylim([-0.005, 0.07]);

%print(gcf, '-depsc', 'ukf_assim_ts.eps');

figure('units','normalized','outerposition',[0 0 1 1])
subplot(311);
plot(repmat(tdays, 1, 2), m_f(:, [4, 5]), 'linewidth', 2);
h = legend('dE', 'dS');
set(h, 'fontsize', 12);
title('Equilibrium changes [UKF]', 'fontsize', 12);

subplot(312);
plot(repmat(tdays, 1, 2), [Ed + m_f(:,4),Ew + m_f(:,4)], 'linewidth', 1.5);
h = legend('Ed+dE', 'Ew+dE');
set(h, 'fontsize', 12);
title('Adjusted equilibria [UKF]', 'fontsize', 12);

subplot(313);
plot(tdays, model_ids(:,2), 'or', 'markerfacecolor', 'red');
title('Active submodel of the moisture model for 10-hr fuel [UKF]', 'fontsize', 12);
ylim([-1, 5]);

%print(gcf, '-depsc', 'ukf_assim_params.eps');
