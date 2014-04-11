


function plot_tsm_loo_errors(sns,errs,aerrs,rmse)

    N = length(sns);
    
    figure;
    ax = subplot(211);
    plot(1:N,errs,'bo-',1:N,aerrs,'ro-',[1;N],[0;0],'k-');
    set(ax,'XTick',[]);
    ylabel('Error [g/100g]');
    
    ax = subplot(212);
    plot(1:N,rmse,'go-');
    ylabel('RMSE');
    set(ax,'XTick',1:N);
    set(ax,'XTickLabel',sns);
    
    
    