
year = '2013';
t_fcast = 48;

[errs1,aerrs1] = forecast_ability_station_data(year,2,t_fcast,1);
[errs2,aerrs2] = forecast_ability_station_data(year,7,t_fcast,1);
[errs3,aerrs3] = forecast_ability_station_data(year,13,t_fcast,1);
[errs4,aerrs4] = forecast_ability_station_data(year,25,t_fcast,1);

[errs5,aerrs5] = forecast_ability_station_data(year,6,t_fcast,2);
[errs6,aerrs6] = forecast_ability_station_data(year,12,t_fcast,2);
[errs7,aerrs7] = forecast_ability_station_data(year,24,t_fcast,2);

[errs9,aerrs9] = forecast_ability_station_data(year,8,t_fcast,4);
[errs10,aerrs10] = forecast_ability_station_data(year,12,t_fcast,4);
[errs11,aerrs11] = forecast_ability_station_data(year,24,t_fcast,4);

[errs12,aerrs12] = forecast_ability_station_data(year,24,t_fcast,12);
[errs13,aerrs13] = forecast_ability_station_data(year,48,t_fcast,12);


% store the results as a regular matrix
steps = [1,1,1,1,2,2,2,4,4,4,12,12];
t_inits = [2,7,13,25,6,12,24,8,12,24,24,48];

save(['station_data_fcast_',year,'.mat'], 'year', 't_fcast', 'steps', 't_inits', ...
     'errs1', 'aerrs1', 'errs2', 'aerrs2', 'errs3', 'aerrs3', 'errs4', 'aerrs4',  ...
     'errs5', 'aerrs5', 'errs6', 'aerrs6', 'errs7', 'aerrs7', ...
     'errs9', 'aerrs9', 'errs10', 'aerrs10', 'errs11', 'aerrs11', 'errs12', 'aerrs12',  ...
     'errs13', 'aerrs13');
     