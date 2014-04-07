
%
% Compute forecast errors for different starting times, where observations
% are assimilated for time t_init (with observations 'step' apart) and then
% the forecast is tested for t_fcast.
%
% synopsis: [err_ekf,err_ukf,err_nf] = estimate_forecast_errors_station(sd,froms,t_init,t_fcast,step)
%

function [err_ekf,err_ukf,err_nf] = estimate_forecast_errors_station(sd,froms,t_init,t_fcast,step)
    N = length(froms);
    M = t_init+t_fcast;

    err_ekf = zeros(N,M);
    err_ukf = zeros(N,M);
    err_nf = zeros(N,M);
    ndx = 1;

    for i=1:N
        rng = froms(i):froms(i)+t_init+t_fcast;
        dts = diff(sd.tdays(rng)) * 24 - 1;
        if(all(abs(dts) < 1/12) && ~any(isnan(sd.fm10(rng))))
            [~,fm10,ekfs,ukfs,nfs] = compare_forecast_vs_kf(sd,froms(i),t_init,t_fcast,step);
            err_ekf(ndx,:) = ekfs - fm10;
            err_ukf(ndx,:) = ukfs - fm10;
            err_nf(ndx,:) = nfs - fm10;
            ndx = ndx + 1;
        end
    end
    
    err_ekf = err_ekf(1:ndx-1,:);
    err_ukf = err_ukf(1:ndx-1,:);
    err_nf = err_nf(1:ndx-1,:);
