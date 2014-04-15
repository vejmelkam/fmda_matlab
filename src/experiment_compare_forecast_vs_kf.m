
year = '2013';
t_fcast = 48;

[errs1,aerrs1] = forecast_ability_station_data(year,2,t_fcast,1);
drawnow;
[errs2,aerrs2] = forecast_ability_station_data(year,7,t_fcast,1);
drawnow;
[errs3,aerrs3] = forecast_ability_station_data(year,13,t_fcast,1);
drawnow;
[errs4,aerrs4] = forecast_ability_station_data(year,25,t_fcast,1);
drawnow;

[errs5,aerrs5] = forecast_ability_station_data(year,6,t_fcast,2);
drawnow;
[errs6,aerrs6] = forecast_ability_station_data(year,12,t_fcast,2);
drawnow;
[errs7,aerrs7] = forecast_ability_station_data(year,24,t_fcast,2);
drawnow;

[errs9,aerrs9] = forecast_ability_station_data(year,8,t_fcast,4);
drawnow;
[errs10,aerrs10] = forecast_ability_station_data(year,12,t_fcast,4);
drawnow;
[errs11,aerrs11] = forecast_ability_station_data(year,24,t_fcast,4);
drawnow;

[errs12,aerrs12] = forecast_ability_station_data(year,24,t_fcast,12);
drawnow;
[errs13,aerrs13] = forecast_ability_station_data(year,48,t_fcast,12);
drawnow;


% store the results as a regular matrix
steps = [1,1,1,1,2,2,2,4,4,4,12,12];
t_inits = [2,7,13,25,6,12,24,8,12,24,24,48];

save(['station_data_fcast_',year,'.mat'], 'year', 't_fcast', 'steps', 't_inits', ...
     'errs1', 'aerrs1', 'errs2', 'aerrs2', 'errs3', 'aerrs3', 'errs4', 'aerrs4',  ...
     'errs5', 'aerrs5', 'errs6', 'aerrs6', 'errs7', 'aerrs7', ...
     'errs9', 'aerrs9', 'errs10', 'aerrs10', 'errs11', 'aerrs11', 'errs12', 'aerrs12',  ...
     'errs13', 'aerrs13');
     
errs_arr = { errs1, errs2, errs3, errs4, errs5, errs6, ...
              errs7, errs9, errs10, errs11, errs12, errs13 };
 
aerrs_arr = { aerrs1, aerrs2, aerrs3, aerrs4, aerrs5, aerrs6, ...
              aerrs7, aerrs9, aerrs10, aerrs11, aerrs12, aerrs13 };

% figure;
% for i=[4,7,9,11]
%     e = errs_arr{i};
%     N = size(e,1);
%     e_ekf = e(:,1); e_ukf = e(:,2);
%     plot(48-N+2:48,e_ekf(2:end),'b-');
%     hold on;
%     plot(48-N+2:48,e_ukf(2:end),'r-');
% end
% hold off;
% title('Average errors');
% 
% figure;
% for i=[4,7,9,11]
%     e = aerrs_arr{i};
%     N = size(e,1);
%     e_ekf = e(:,1); e_ukf = e(:,2);
%     plot(48-N+2:48,e_ekf(2:end),'b-');
%     hold on;
%     plot(48-N+2:48,e_ukf(2:end),'r-');
% end
% hold off;
% title('Absolute errors');
 