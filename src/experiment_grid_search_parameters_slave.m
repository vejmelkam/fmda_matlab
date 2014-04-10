    


function [sse] = experiment_grid_search_parameters_slave(from,to,year)



    yr = num2str(year);
    Lst = dir(['../data/*fm10*',yr,'*']);

    % load data
    sds = cell(length(Lst),1);
    ndx = 1;
    for i=1:length(Lst)
        data = load_station_data(Lst(i).name(1:5),yr,[100,270]);
        if(length(data.tdays) > 10)
            sds{ndx} = data;
            ndx = ndx + 1;
        end
    end
    Nst = ndx - 1;
    fprintf('Loaded %d stations.\n', Nst);
    sds = sds(1:Nst);
    
    data = load('grid_search_items');
    items = data.items;
    items = items(from:to,:);
    fprintf('Beginning parameter search for %d items.\n', size(items,1));
    sse = zeros(size(items,1),1);

    for j=1:size(items,1)
        S = items(j,1);
        rk = items(j,2);
        r0 = items(j,3);
        dE = items(j,4);
        Trk = items(j,5);
        fprintf('S=%g rk=%g r0=%g dE=%g Trk=%g', S,rk,r0,dE,Trk);
        for i=1:Nst
            fprintf('.');
            Ed = sds{i}.ed;
            Ew = sds{i}.ew;
            fm10 = sds{i}.fm10;
            rain = sds{i}.rain;
            tdays = sds{i}.tdays - datenum(year,1,1);

            % time in hours (integration step for moisture is 1h)
            N = length(tdays);

            % parameters of the simulation
            k=3;
            m = zeros(k+2,1);
            m(1:k) = 0.5 * (Ed(1)+Ew(1));
            m(k+1) = dE;

            m_t = zeros(N, 5);
            m_t(1, :) = m';
            for t=2:N
                dt = (tdays(t) - tdays(t-1)) * 86400;
                Ed2 = (Ed(t)+Ed(t-1))/2;
                Ew2 = (Ew(t)+Ew(t-1))/2;

                m_new = moisture_model_ext2([1,10,100]',Ed2,Ew2,m,rain(t),dt,1e10,S,rk,r0,Trk);
                m_t(t, :) = m_new;
                m = m_new;
            end

            % compare values
            ndx = isfinite(fm10);
            sse(j) = sse(j) + sum((fm10(ndx)-m_t(ndx,2)).^2);

        end
        fprintf(' [%g]\n', sse(j));
    end
end
