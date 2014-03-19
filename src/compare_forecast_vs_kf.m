
%
%  Compares the predictive ability resulting from data assimilation as
%  follows.  A moisture model starts from equilibrium and is run for
%  't_init' steps in these three configurations:
%
%    - no assimilation
%    - EKF assimilation
%    - UKF assimilation
%
%  After 't_init' the moisture model is run simply in forecast mode for
%  't_fcast' time steps and the difference between future (unseen during
%  t_init phase) observations and the model forecast is evaluated.
%
%

function [tm,fm10s,m_ekfs,m_ukfs,m_nfs] = compare_forecast_vs_kf(sd,from,t_init,t_fcast)

    % pre-compute equilibrium moisture from station data
    [Ed,Ew] = equilibrium_moisture2(sd.relh,sd.t2);
    tdays = sd.tdays - sd.tdays(1);
    tm = tdays(from+1:from+t_init+t_fcast);

    % number of fuels and extended model size
    k = 3;
    M = k+2;
    Npts = 2*M+1;

    % parameters of the simulation
    m_ext = zeros(M,1);
    m_ext(1:k) = 0.5 * (Ed(from+1)+Ew(from+1));

    % initialize covariance
    P_ukf = eye(M) * 0.01;   % error covariance of the initial guess UKF
    P_ekf = eye(M) * 0.01;   % error covariance of the initial guess EKF

    % Kalman filter Q (model error covariance)
    Qphr = zeros(M);
    Qphr(1:k, 1:k) = diag([0.0005,0.00005,0.00001]);
    Qphr(k+1,k+1) = 0.0001;
    Qphr(k+2,k+2) = 0.0001;

    % the observation operator is a n_k x Ndim matrix with I_(n_k) on the left
    H = zeros(1,k+2);
    H(2) = 1;

    % storage space for results (with filtering)
    m_ekf = m_ext;
    m_ukf = m_ext;
    m_nf  = m_ext;

    % allocate space for results
    fm10s = zeros(t_init+t_fcast,1);
    m_ekfs = zeros(t_init+t_fcast,1);
    m_ukfs = zeros(t_init+t_fcast,1);
    m_nfs = zeros(t_init+t_fcast,1);
    
    m_ekfs(1) = m_ekf(2);
    m_ukfs(1) = m_ukf(2);
    m_nfs(1) = m_nf(2);
    fm10s(1) = sd.fm10(from+1);
    
    for i=2:t_init
        
        % compute real world inputs
        dt = (tdays(from+i) - tdays(from+i-1)) * 86400;
        Ed2 = (Ed(from+i)+Ed(from+i-1))/2;
        Ew2 = (Ew(from+i)+Ew(from+i-1))/2;

        % compute Jacobian for EKF
        Jm = moisture_tangent_model_ext([1,10,100]',Ew2,Ed2,m_ekf,sd.rain(from+i),dt,1e10);
        
        % run forecasts for EKF and for (no filtering) NF
        m_ekf = moisture_model_ext([1,10,100]',Ed2,Ew2,m_ekf,sd.rain(from+i),dt,1e10);
        m_nf  = moisture_model_ext([1,10,100]',Ed2,Ew2,m_nf,sd.rain(from+i),dt,1e10);

        % run forecast for UKF
        [m_sigma, w] = ukf_select_sigma_points(m_ukf,P_ukf,0);
        m_sigma_i = zeros(M, Npts);
        for n=1:Npts-1
            m_sigma_i(:,n) = moisture_model_ext([1,10,100]',Ed2,Ew2,m_sigma(:,n),sd.rain(from+i),dt,1e10);
        end
        m_sigma_i(:,Npts) = moisture_model_ext([1,10,100]',Ed2,Ew2,m_ukf,sd.rain(from+i),dt,1e10);

        % compute the prediction mean x_mean(i|i-1) for UKF
        m_pred = sum(m_sigma_i * diag(w), 2);

        % propagate  covariance through Jacobian (EKF)
        P_ekf = Jm*P_ekf*Jm' + Qphr*dt/3600;
        
        % propagate covariance (UKF)
        sqrtP = (m_sigma_i - repmat(m_pred, 1, Npts)) * diag(sqrt(w));
        P_ukf = Qphr*dt/3600 + sqrtP * sqrtP';        

        % retrieve observation
        m_measured = sd.fm10(from+i);
        R = sd.fm10_var(from+i);        
        
        % do assimilation step for EKF
        S = (H*P_ekf*H' + R);
        K_ekf = P_ekf * H' / S;
        m_ekf = m_ekf + K_ekf*(m_measured - H*m_ekf);
        P_ekf = P_ekf - K_ekf*S*K_ekf';
        
        % do assimilation step for UKF
        y_pred_i = H * m_sigma_i;
        y_pred = sum(y_pred_i * diag(w), 2);
        sqrtS = (y_pred_i - repmat(y_pred, 1, Npts)) * diag(sqrt(w));
        S = sqrtS * sqrtS' + R; 
        Cxy = (m_sigma_i - repmat(m_pred, 1, Npts)) * diag(w) * (y_pred_i - repmat(y_pred, 1, Npts))';
        K_ukf = Cxy / S;
        m_ukf = m_pred + K_ukf*(m_measured - y_pred);
        P_ukf = P_ukf - K_ukf*S*K_ukf';
        
        m_ekfs(i) = m_ekf(2);
        m_ukfs(i) = m_ukf(2);
        m_nfs(i) = m_nf(2);
        fm10s(i) = sd.fm10(from+i);
        
    end
        
    for i=t_init+1:t_init+t_fcast
        
        % compute real world inputs
        dt = (tdays(from+i) - tdays(from+i-1)) * 86400;        
        Ed2 = (Ed(from+i)+Ed(from+i-1))/2;
        Ew2 = (Ew(from+i)+Ew(from+i-1))/2;
        
        m_ekf = moisture_model_ext([1,10,100]',Ed2,Ew2,m_ekf,sd.rain(from+i),dt,1e10);
        m_ukf = moisture_model_ext([1,10,100]',Ed2,Ew2,m_ukf,sd.rain(from+i),dt,1e10);
        m_nf = moisture_model_ext([1,10,100]',Ed2,Ew2,m_nf,sd.rain(from+i),dt,1e10);
        
        m_ekfs(i) = m_ekf(2);
        m_ukfs(i) = m_ukf(2);
        m_nfs(i) = m_nf(2);
        fm10s(i) = sd.fm10(from+i);
        
    end
    
        
        
        
    