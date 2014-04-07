

%
%  This script runs the extended moisture model for one grid point with the
%  Unscented Kalman Filter.
%

function run_moisture_model_ext_with_ukf(station,yr)

    % Load station data & compute equilibria
    sd = load_station_data(station,yr);
    tdays = sd.tdays - sd.tdays(1);
    Ed = sd.ed;
    Ew = sd.ew;

    % time in hours (integration step for moisture is 1h)
    N = length(sd.tdays);
    do_assim = false(N,1);
    do_assim(1:6:N) = true;

    % number of fuels and extended model size
    k = 3;
    Tk = [1,10,100]';
    M = k+2;

    % parameters of the simulation
    m_ext = zeros(M,1);
    m_ext(1:k) = 0.5 * (Ed(1)+Ew(1));

    % initialize covariance
    P = eye(M) * 0.01;   % error covariance of the initial guess

    % Kalman filter Q (model error covariance)
    Qphr = zeros(M);
    Qphr(1:k, 1:k) = diag([0.0005,0.00005,0.00001]);
    Qphr(k+1,k+1) = 0.0001;
    Qphr(k+2,k+2) = 0.0001;

    % the observation operator is a n_k x Ndim matrix with I_(n_k) on the left
    H = zeros(1,k+2);
    H(2) = 1;

    % storage space for results (with filtering)
    m_f = zeros(N, M);
    m_f(1, :) = m_ext';

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
    W0 = 1/(Npts+1);

    % Storage of the sigma points
    sigma_pts = zeros(N, M, Npts);

    % predict & update loop
    for i = 2:N

        % compute the integration time step
        dt = (sd.tdays(i) - sd.tdays(i-1)) * 86400;
        Ed2 = (Ed(i)+Ed(i-1))/2;
        Ew2 = (Ew(i)+Ew(i-1))/2;

        % draw the sigma points
        [m_sigma, w] = ukf_select_sigma_points(m_f(i-1,:)', P, W0);
        sigma_pts(i, :, :) = m_sigma;

        % compute & store results for system without KF
        m_n(i, :) = moisture_model_ext(Tk,Ed2,Ew2,m_n(i-1,:)',sd.rain(i),dt,1e10);

        % compute & store results for corrected system without KF
        m_c(i, :) = moisture_model_ext(Tk,Ed2,Ew2,m_c(i-1,:)',sd.rain(i),dt,1e10);

        % UKF forecast
        [mf,sqrtP,Pf] = ukf_forecast(Tk,Ed2,Ew2,m_f(i-1,:)',sd.rain(i),dt,1e10,P,Qphr);

        % KALMAN UPDATE STEP (if measurement available) 
        if(do_assim(i) && ~isnan(sd.fm10(i)))

            [m_f(i,:),P,K,S] = ukf_update(mf,sqrtP,Pf,H,sd.fm10(i),sd.fm10_var(i));

            trS(i) = prod(eig(S));
            sK(i,:) = K;
            trP(i) = abs(prod(eig(P)));

        else

            % if no observation is available, store the predicted value
            m_f(i,:) = mf;

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
    fprintf('MSE for UKF: %g  MAPE for UKF: %g\n', norm(sd.fm10(ndx) - m_f(ndx,2),2)/N,norm(sd.fm10(ndx) - m_f(ndx,2),1)/N);
    
    % store results
    savefig(f, [station,'_',yr,'_ukf.fig']);
    save([station,'_',yr,'_ukf.mat'], 'sd', 'm_f','m_n');

end
