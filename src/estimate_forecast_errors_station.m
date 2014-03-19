


function [err_ekf,err_ukf,err_nf] = estimate_forecast_errors_station(sd,froms,t_init,t_fcast)
    N = length(froms);
    M = t_init+t_fcast;

    err_ekf = zeros(N,M);
    err_ukf = zeros(N,M);
    err_nf = zeros(N,M);
    ndx = 1;

    for i=1:N
        dts = diff(sd.tdays(froms(i):froms(i)+t_init+t_fcast)) * 24 - 1;
        if(all(abs(dts) < 1/12))
            [~,fm10,ekfs,ukfs,nfs] = compare_forecast_vs_kf(sd,froms(i),t_init,t_fcast);
            err_ekf(ndx,:) = ekfs - fm10;
            err_ukf(ndx,:) = ukfs - fm10;
            err_nf(ndx,:) = nfs - fm10;
            ndx = ndx + 1;
        end
    end
    
    err_ekf = err_ekf(1:ndx-1,:);
    err_ukf = err_ukf(1:ndx-1,:);
    err_nf = err_nf(1:ndx-1,:);
