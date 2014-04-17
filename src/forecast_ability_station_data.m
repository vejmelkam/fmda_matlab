
%
%  An experiment where the forecasting skill of UKF/EKF/no filter
%  initialized model are compared.  This test is aggregated over all
%  stations.
%
%
%  synopsis: [err_ekf,err_ukf,err_nf,stations] = experiment_forecast_ability_station_data(year,t_init,t_fcast,step)
%
%

function [err_ekf,err_ukf,err_nf,stations] = forecast_ability_station_data(year,t_init,t_fcast,step)

    S = dir(['../data/*fm10*',year,'*']);
    N = length(S);
    
    stavg = zeros(t_init+t_fcast,N,3);
    stabs = zeros(t_init+t_fcast,N,3);

    for i=1:N

        station = S(i).name(1:5);
        fprintf('%s ', station);
        if(rem(i,10) == 0), fprintf('\n'); end
        sd = load_station_data(station, year);

        if(i==1)
            [err_ekf,err_ukf,err_nf] = estimate_forecast_errors_station(sd,1:9:length(sd.tdays)-t_fcast-t_init,t_init,t_fcast,step);
            stations = ones(size(err_ekf,1),1);
            stavg(:,i,:) = [mean(err_ekf);mean(err_ukf);mean(err_nf)]';
            stabs(:,i,:) = [mean(abs(err_ekf));mean(abs(err_ukf));mean(abs(err_nf))]';
        else
            [err_ekf_i,err_ukf_i,err_nf_i] = estimate_forecast_errors_station(sd,1:9:length(sd.tdays)-t_fcast-t_init,t_init,t_fcast,step);
            err_ekf = [err_ekf;err_ekf_i]; %#ok<AGROW>
            err_ukf = [err_ukf;err_ukf_i]; %#ok<AGROW>
            err_nf  = [err_nf;err_nf_i]; %#ok<AGROW>
            stations = [stations;i*ones(size(err_ekf_i,1),1)]; %#ok<AGROW>
            stavg(:,i,:) = [mean(err_ekf_i);mean(err_ukf_i);mean(err_nf_i)]';
            stabs(:,i,:) = [mean(abs(err_ekf_i));mean(abs(err_ukf_i));mean(abs(err_nf_i))]';
        end

    end
    fprintf('\n');
    
    % transpose all to have time in rows
    err_ekf = err_ekf';
    err_ukf = err_ukf';
    err_nf = err_nf';

    figure;
    subplot(241);
    plot(stavg(:,:,1));
    title(sprintf('EKF [t_init=%d,t_fcast=%d]', t_init,t_fcast));
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(245);
    plot(stabs(:,:,1));
    title('EKF');
    xlabel('Time [hours]');
    ylabel('Abs. error [-]');

    subplot(242);
    plot(stavg(:,:,2));
    title('UKF');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(246);
    plot(stabs(:,:,2));
    title('UKF');
    xlabel('Time [hours]');
    ylabel('Abs. Error [-]');

    subplot(243);
    plot(stavg(:,:,3));
    title('None');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(247);
    plot(stabs(:,:,3));
    title('None');
    xlabel('Time [hours]');
    ylabel('Abs. Error [-]');

    subplot(244);
    plot([nanmean(err_ekf,2),nanmean(err_ukf,2),nanmean(err_nf,2)]);
%    plot([nanmean(stavg(:,:,1),2),nanmean(stavg(:,:,2),2),nanmean(stavg(:,:,3),2)]);
    legend('ekf','ukf','none');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(248);
    plot([nanmean(abs(err_ekf),2),nanmean(abs(err_ukf),2),nanmean(abs(err_nf),2)]);
%    plot([nanmean(stabs(:,:,1),2),nanmean(stabs(:,:,2),2),nanmean(stabs(:,:,3),2)]);
    legend('ekf','ukf','none');
    xlabel('Time [hours]');
    ylabel('Error [-]');
    
    drawnow;
    
end

