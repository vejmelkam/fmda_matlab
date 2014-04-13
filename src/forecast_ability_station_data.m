
%
%  An experiment where the forecasting skill of UKF/EKF/no filter
%  initialized model are compared.  This test is aggregated over all
%  stations.
%
%
%  synopsis: [errs,aerrs] = experiment_forecast_ability_station_data(year,t_init,t_fcast,step)
%
%

function [errs,aerrs] = forecast_ability_station_data(year,t_init,t_fcast,step)

    S = dir(['../data/*fm10*',year,'*']);
    N = length(S);

    all_ekfs = zeros(N,t_init+t_fcast);
    all_ukfs = zeros(N,t_init+t_fcast);
    all_nfs = zeros(N,t_init+t_fcast);

    all_abs_ekfs = zeros(N,t_init+t_fcast);
    all_abs_ukfs = zeros(N,t_init+t_fcast);
    all_abs_nfs = zeros(N,t_init+t_fcast);

    for i=1:length(S)

        station = S(i).name(1:5);
        fprintf('%s ', station);
        if(rem(i,10) == 0), fprintf('\n'); end
        sd = load_station_data(station, year);

        [err_ekf,err_ukf,err_nf] = estimate_forecast_errors_station(sd,1:7:length(sd.tdays)-t_fcast-t_init,t_init,t_fcast,step);

        all_ekfs(i,:) = mean(err_ekf);
        all_ukfs(i,:) = mean(err_ukf);
        all_nfs(i,:) = mean(err_nf);

        all_abs_ekfs(i,:) = mean(abs(err_ekf));
        all_abs_ukfs(i,:) = mean(abs(err_ukf));
        all_abs_nfs(i,:) = mean(abs(err_nf));

    end
    fprintf('\n');

    figure;
    subplot(241);
    plot(all_ekfs');
    title(sprintf('EKF [t_init=%d,t_fcast=%d]', t_init,t_fcast));
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(245);
    plot(all_abs_ekfs');
    title('EKF');
    xlabel('Time [hours]');
    ylabel('Abs. error [-]');

    subplot(242);
    plot(all_ukfs');
    title('UKF');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(246);
    plot(all_abs_ukfs');
    title('UKF');
    xlabel('Time [hours]');
    ylabel('Abs. Error [-]');

    subplot(243);
    plot(all_nfs');
    title('None');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(247);
    plot(all_abs_nfs');
    title('None');
    xlabel('Time [hours]');
    ylabel('Abs. Error [-]');

    subplot(244);
    plot([nanmean(all_ekfs)',nanmean(all_ukfs)',nanmean(all_nfs)']);
    legend('ekf','ukf','none');
    xlabel('Time [hours]');
    ylabel('Error [-]');

    subplot(248);
    plot([nanmean(all_abs_ekfs)',nanmean(all_abs_ukfs)',nanmean(all_abs_nfs)']);
    legend('ekf','ukf','none');
    xlabel('Time [hours]');
    ylabel('Error [-]');
    
    drawnow;

    eamean = nanmean(all_abs_ekfs)';
    uamean = nanmean(all_abs_ukfs)';
    emean = nanmean(all_ekfs)';
    umean = nanmean(all_ukfs)';

    errs = [emean,umean];
    aerrs = [eamean,uamean];

end

