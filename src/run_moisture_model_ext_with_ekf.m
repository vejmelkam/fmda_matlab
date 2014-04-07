
%
%  This script runs the extended moisture model for one grid point with the
%  Extended Kalman filter.
%

function run_moisture_model_ext_with_ekf(station,yr)

% Load station data & compute equilibria
sd = load_station_data(station,yr);
tdays = sd.tdays - sd.tdays(1);
Ed = sd.ed;
Ew = sd.ew;

% time in hours (integration step for moisture is 1h)
N = length(tdays);
do_assim = false(N,1);
do_assim(1:12:N) = true;

% number of fuels and extended model size
Tk = [1,10,100]';
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

% predict & update 0loop
for i=2:N
    
    % compute the integration time step in seconds
    dt = (tdays(i) - tdays(i-1)) * 86400;
    Ed2 = (Ed(i)+Ed(i-1))/2;
    Ew2 = (Ew(i)+Ew(i-1))/2;
    
    % compute & store results for system without Kalman filtering
    m_n(i, :) = moisture_model_ext(Tk,Ed2,Ew2,m_n(i-1,:)',sd.rain(i),dt,1e10);
    
    % compute & store results for system without Kalman filtering but with
    % initial dE & dS correction
    m_c(i, :) = moisture_model_ext(Tk,Ed2,Ew2,m_c(i-1,:)',sd.rain(i),dt,1e10);

    % KALMAN PREDICT STEP
    % estimate new moisture mean based on last best guess (m)
    [m_pred,P,J,model_ids(i,:)] = ekf_forecast(Tk,Ed2,Ew2,m_f(i-1,:)',sd.rain(i),dt,1e10,P,Qphr);
    sP(i, :, :) = P;
    trP(i) = prod(eig(P(1:k,1:k)));
    trJ(i) = prod(eig(J));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if(do_assim(i) && ~isnan(sd.fm10(i)))
        
        [m_f(i,:),P,K,S] = ekf_update(m_pred,P,H,sd.fm10(i),sd.fm10_var(i));
        trS(i) = prod(eig(S));
        sK(i,:) = K;
    
    else
        
        % if no observation is available, store the predicted value
        m_f(i,:) = m_pred;
        
    end
        
end

    set(0,'DefaultAxesLooseInset',[0,0,0,0])

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(411);
    plot(tdays, m_f(:,2), 'r-', 'linewidth', 1.5);
    hold on;
    plot(tdays, m_n(:,2), 'g-', 'linewidth', 1.5);
    %plot(tdays, m_c(:,2), 'b-', 'linewidth', 1.5);
    %plot(tdays, rain, 'k--', 'linewidth', 1.5);
    plot(tdays(do_assim), sd.fm10(do_assim), 'ko', 'markersize', 4, 'markerfacecolor', 'm');
    %plot(repmat(tdays, 1, 2), [m_f(:,2) - sqrt(sP(:, 2, 2)), m_f(:,2) + sqrt(sP(:, 2, 2))], 'rx');
    h = legend('sys + UKF', 'sys', 'obs', 'orientation', 'horizontal');
    set(h, 'fontsize', 10);
    title(['Evolution of the moisture model [',station,'] [UKF]'], 'fontsize', 12);
    ylim([0, min(1.2,1.1*max([m_f(:,2);m_n(:,2);m_c(:,2)]))]);

    % select time indices corresponding to observation times
    subplot(412);
    % plot(tdays, log10(trP), 'b-', 'linewidth', 2);
    %plot(tdays(do_assim), log10(trS(do_assim)), 'ko', 'markerfacecolor', 'green', 'markersize', 4);
    plot(tdays(do_assim), log10(abs(sK(do_assim,1))), 'ko', 'markerfacecolor', 'red', 'markersize', 4);
    hold on;
    plot(tdays(do_assim), log10(abs(sK(do_assim,2))), 'ko', 'markerfacecolor', 'green', 'markersize', 4);
    plot(tdays(do_assim), log10(abs(sK(do_assim,3))), 'ko', 'markerfacecolor', 'blue', 'markersize', 4);
    hold off;
    h = legend('K(1)', 'K(2)', 'K(3)', 'orientation', 'horizontal');
    set(h, 'fontsize', 10);
    title('Kalman filter: Kalman gains for 1-hr,10-hr,100-hr fuel [UKF]', 'fontsize', 12);
    ylim([-12, 5]);

%     subplot(313);
%     plot(tdays, sP(:, 2, 2), 'r-', 'linewidth', 2);
%     hold on
%     plot(tdays, sP(:, 2, 4), 'g-', 'linewidth', 2);
%     plot(tdays, sP(:, 2, 5), 'b-', 'linewidth', 2);
%     plot(tdays, sP(:, 1, 1), 'k-', 'linewidth', 2);
%     plot(tdays, sP(:, 3, 3), 'c--', 'linewidth', 2);
%     hold off
%     h = legend('var(m(2))', 'cov(m(2),dE)', 'cov(m(2),dS)', 'var(m(1))', 'var(m(3))', 'orientation', 'horizontal');
%     set(h, 'fontsize', 12);
%     title('Covariance between moisture and system parameters [UKF]', 'fontsize', 12);

    %print(gcf, '-depsc', 'ukf_assim_ts.eps');

%    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(413);
    plot(repmat(tdays, 1, 2), m_f(:, [4, 5]), 'linewidth', 2);
    h = legend('dE', 'dS');
    set(h, 'fontsize', 12);
    title('Equilibrium changes [UKF]', 'fontsize', 12);

    subplot(414);
    plot(repmat(tdays, 1, 2), [Ed + m_f(:,4),Ew + m_f(:,4)], 'linewidth', 1.5);
    h = legend('Ed+dE', 'Ew+dE');
    set(h, 'fontsize', 12);
    title('Adjusted equilibria [UKF]', 'fontsize', 12);
%     subplot(413);
%     plot(tdays, model_ids(:,2), 'or', 'markerfacecolor', 'red');
%     title('Active submodel of the moisture model for 10-hr fuel [UKF]', 'fontsize', 12);
%     ylim([-1, 5]);
% 
%     subplot(414);
%     plot(tdays(do_assim), [sd.fm10(do_assim),sd.relh(do_assim),sd.rain(do_assim)]);
%     title('Input data fuel model [UKF]', 'fontsize', 12);
%     legend('fm10','relh','rain');
% 
    ndx = isfinite(sd.fm10);
    fprintf('MSE for EKF: %g  MAPE for EKF: %g\n', norm(sd.fm10(ndx) - m_f(ndx,2),2)/N,norm(sd.fm10(ndx) - m_f(ndx,2),1)/N);
    
    % store results
    savefig(f, [station,'_',yr,'_ekf.fig']);
    save([station,'_',yr,'_ekf.mat'], 'sd', 'm_f','m_n');

end
