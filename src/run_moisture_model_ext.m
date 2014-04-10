

%
%  This script runs the extended moisture model for one spatial point
%

function run_moisture_model_ext(station,yr)
    % Load station data
    sd = load_station_data(station,yr,[100,270]);
    Ed = sd.ed;
    Ew = sd.ew;
    sd.tdays = sd.tdays - datenum(2012,1,1);
    
    if(isempty(Ed))
        return;
    end

    % time in hours (integration step for moisture is 1h)
    N = length(sd.tdays);

    % parameters of the simulation
    k=3;
    m = zeros(k+2,1);
    m(1:k) = 0.5 * (Ed(1)+Ew(1));

    m_t = zeros(N, 5);
    m_t(1,:) = m';
    m_p = zeros(N, 5);
    m_p(1,:) = m';
    m_p(1,4) = -0.04;
    for i=2:N
        dt = (sd.tdays(i) - sd.tdays(i-1)) * 86400;
        Ed2 = (Ed(i)+Ed(i-1))/2;
        Ew2 = (Ew(i)+Ew(i-1))/2;

        m_new = moisture_model_ext([1,10,100]',Ed2,Ew2,m,sd.rain(i),dt,1e10);
        m_t(i, :) = m_new;
        m = m_new;
        
        m_p(i,:) = moisture_model_ext2([1,10,100]',Ed2,Ew2,m_p(i-1,:)',sd.rain(i),dt,1e10,0.6,2,0.04,9);
    end

    f = figure();
    subplot(4,1,1);
    plot(repmat(sd.tdays, 1, 3), [m_t(:,2),m_p(:,2),sd.fm10], 'linewidth', 1.25);
    ylabel('10-hr fuel moisture [g/100g]');
    title(['Model vs. observations in ',yr,' station: ', station], 'fontsize', 14);
    legend('model orig.', 'model fitted', 'observations');
    axis tight
    subplot(4,1,2);
    plot(repmat(sd.tdays, 1, 2), [Ew,Ed], 'linewidth', 1.25);
    ylabel('model equilibria [g/100g]');
    legend('wetting', 'drying');
    axis tight
    subplot(4,1,3);
    plot(sd.tdays,sd.rain,'linewidth',1.25);
    ylabel('rain intensity [mm/h]');
    xlabel('days from start of year');
    axis tight
    subplot(4,4,13);
    scatter(sd.fm10,m_t(:,2),3,'markerfacecolor','b');
    hold on;
    plot([0, max([m_t(:,2);sd.fm10])], [0,max([m_t(:,2);sd.fm10])], 'k-');
    hold off;
    axis square tight
    xlabel('fm10 obs [g/100g]');
    ylabel('fm10 model orig [g/100g]');
    subplot(4,4,14);
    scatter(sd.fm10,m_p(:,2),3,'markerfacecolor','b');
    hold on;
    plot([0, max([m_t(:,2);sd.fm10])], [0,max([m_t(:,2);sd.fm10])], 'k-');
    hold off;
    axis square tight
    xlabel('fm10 obs [g/100g]');
    ylabel('fm10 model fitted [g/100g]');
    subplot(4,4,15);
    scatter(sd.fm10,sd.t2,3,'markerfacecolor','b');
    xlabel('fm10 obs [g/100g]');
    ylabel('t2 obs [K]');
    axis square tight
    subplot(4,4,16);
    scatter(sd.fm10,sd.relh,3,'markerfacecolor','b');
    xlabel('fm10 obs [g/100g]');
    ylabel('relh obs [-]');
    axis square tight  

    savefig(f, [station,'-',yr,'_model_vs_obs.fig']);
    save([station,'-',yr,'_model_vs_obs.mat'],'sd','m_t','m_p');
end

