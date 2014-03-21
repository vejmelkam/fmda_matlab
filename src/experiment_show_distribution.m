

%
%
%  In this experiment the histogram of observations of a certain stations.
%
%

function experiment_show_distribution(station_code,year,R)

    if(nargin < 3)
        R = [];
    end

    sd = load_station_data(station_code,year,R);
    
    Nt = length(sd.tdays);
    have_fm10 = isfinite(sd.fm10);
    Nf = sum(have_fm10);
    fprintf('No. of observations: %d (fm10: %d)\n', ...
            Nt,Nf);
        
    figure;
    
    subplot(2,3,1);
    hist(sd.fm10, ceil(sqrt(Nf)));
    title('fm10 observations');
    xlabel('fraction [-]');
    
    subplot(2,3,2);
    hist(sd.relh, ceil(sqrt(Nt)));
    title('rel. humidity observations');
    xlabel('fraction [-]');

    subplot(2,3,3);
    hist(sd.rain(sd.rain > 0), ceil(sqrt(Nt)));
    xlabel('[mm/h]');
    title(sprintf('rain observations [non-zero only], zeros: %d', sum(sd.rain==0)));

    subplot(2,3,4);
    hist(sd.t2, ceil(sqrt(Nt)));
    title('T2 observations');
    xlabel('K');
    
    subplot(2,3,5);
    scatter(sd.relh(have_fm10),sd.fm10(have_fm10));
    title('Scatter: fm10 vs. relh');
    xlabel('relh [-]');
    ylabel('fm10 [-]');
    
    subplot(2,3,6);
    scatter(sd.t2(have_fm10),sd.fm10(have_fm10));
    title('Scatter: fm10 vs. t2');
    xlabel('T2 [K]');
    ylabel('fm10 [-]');
    
    fm10_mx = max(sd.fm10(have_fm10));
    fprintf('%s: max fm10: %g, count in dataset %d, zero fm10 %d\n', sd.stid, fm10_mx,sum(sd.fm10==fm10_mx),sum(sd.fm10==0));
    

    