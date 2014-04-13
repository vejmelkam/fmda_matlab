
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

function [tm,fm10s,m_ekfs,m_ukfs,m_nfs] = compare_forecast_vs_kf(sd,from,t_init,t_fcast,step)

    % pre-compute equilibrium moisture from station data
    [Ed,Ew] = equilibrium_moisture2(sd.relh,sd.t2);
    tdays = sd.tdays - sd.tdays(1);
    tm = tdays(from+1:from+t_init+t_fcast);

    % number of fuels and extended model size
    Tk = [1,10,100]';
    k = 3;
    M = k+2;

    % parameters of the simulation
    m_ext = zeros(M,1);
    m_ext(1:k) = 0.5 * (Ed(from+1)+Ew(from+1));
    m_ext(k+1) = -0.04; % Colorado 2012 best fit

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
    
    kappa = 0;
    
    for i=2:t_init
        
        % compute real world inputs
        dt = (tdays(from+i) - tdays(from+i-1)) * 86400;
        Ed2 = (Ed(from+i)+Ed(from+i-1))/2;
        Ew2 = (Ew(from+i)+Ew(from+i-1))/2;

        % retrieve observation
        d = sd.fm10(from+i);
        R = sd.fm10_var(from+i);
        
        % run forecasts (uses Colorado 2012 best fits)
        m_nf  = moisture_model_ext2(Tk,Ed2,Ew2,m_nf,sd.rain(from+i),dt,1e10,0.6,2,0.08,7);
        [m_ekf,P_ekf] = ekf_forecast2(Tk,Ed2,Ew2,m_ekf,sd.rain(from+i),dt,1e10,P_ekf,Qphr,0.6,2,0.08,7);
        f = @(x,w) moisture_model_ext2(Tk,Ed2,Ew2,x,sd.rain(from+i),dt,1e10,0.6,2,0.08,7) + w;
        [m_ukf,sqrtP,sigma_f] = ukf_forecast_general(m_ukf,f,P_ukf,Qphr*dt/3600,1,kappa);
        
        if(rem(i,step)==0)
            [m_ekf,P_ekf] = ekf_update(m_ekf,P_ekf,H,d,R);
            [m_ukf,P_ukf] = ukf_update(m_ukf,sqrtP,sigma_f,H,d,R,kappa);
        end
        
        m_ekfs(i) = m_ekf(2);
        m_ukfs(i) = m_ukf(2);
        m_nfs(i) = m_nf(2);
        fm10s(i) = d;
        
    end
        
    for i=t_init+1:t_init+t_fcast
        
        % compute real world inputs
        dt = (tdays(from+i) - tdays(from+i-1)) * 86400;        
        Ed2 = (Ed(from+i)+Ed(from+i-1))/2;
        Ew2 = (Ew(from+i)+Ew(from+i-1))/2;
        
        m_ekf = moisture_model_ext2(Tk,Ed2,Ew2,m_ekf,sd.rain(from+i),dt,1e10,0.6,2,0.08,7);
        m_ukf = moisture_model_ext2(Tk,Ed2,Ew2,m_ukf,sd.rain(from+i),dt,1e10,0.6,2,0.08,7);
        m_nf = moisture_model_ext2(Tk,Ed2,Ew2,m_nf,sd.rain(from+i),dt,1e10,0.6,2,0.08,7);
        
        m_ekfs(i) = m_ekf(2);
        m_ukfs(i) = m_ukf(2);
        m_nfs(i) = m_nf(2);
        fm10s(i) = sd.fm10(from+i);
        
    end
    
        
        
        
    