

%
%  This script runs the extended moisture model for one spatial point
%

% Load ESPC2 data
station = 'ESPC2';
yr = '2012';
[tdays,t2,relh,fm10,rain,accp] = load_station_data(station,yr);

[Ed,Ew] = equilibrium_moisture2(relh,t2);

% time in hours (integration step for moisture is 1h)
N = length(tdays);

% parameters of the simulation
k=3;
m = zeros(k+2,1);
m(1:k) = 0.5 * (Ed(1)+Ew(1));

m_t = zeros(N, 5);
m_t(1, :) = m';
for i=2:N
    dt = (tdays(i) - tdays(i-1)) * 86400;
    m_new = moisture_model_ext([1,10,100]',Ed(i),Ew(i),m,rain(i),dt,1e10);
    m_t(i, :) = m_new;
    m = m_new;
end

subplot(3,1,1);
plot(repmat(tdays, 1, 2), [m_t(:, 2),fm10], 'linewidth', 1.5);
ylabel('fuel moisture [-]');
legend('model', 'observations');
subplot(3,1,2);
plot(repmat(tdays, 1, 2), [Ew,Ed], 'linewidth', 1.5);
ylabel('equilibria [-]');
legend('wetting', 'drying');
subplot(3,1,3);
plot(tdays,rain,'linewidth',1.5);
ylabel('rain');
xlabel(station);

