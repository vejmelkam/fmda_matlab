    


function [sse,Ss,dEs,rks,r0s,Trks] = experiment_grid_search_parameters(year)

    yr = num2str(year);
    Lst = dir(['../data/*fm10*',yr,'*']);
    
    done = false;

    % initial grid search
    dEs = -0.05:0.01:-0.03;
    Ss = 0.2:0.2:0.6;
    rks = 7:1:9;
    r0s = 0.06:0.02:0.1;
    Trks = 9:1:11;

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
    sse = zeros(length(Ss),length(rks),length(r0s),length(dEs),length(Trks));

    while(~done)
        
        fprintf('Limits: S in [%g,%g], dE in [%g,%g], rk in [%d,%d], r0 in [%g,%g], Trk in [%d,%d]\n', ...
                Ss(1),Ss(end),dEs(1),dEs(end),rks(1),rks(end),r0s(1),r0s(end),Trks(1),Trks(end));
        for sndx = 1:length(Ss)
            S = Ss(sndx);
            for rkndx = 1:length(rks)
                rk = rks(rkndx);
                for r0ndx = 1:length(r0s)
                    r0 = r0s(r0ndx);
                    for dendx = 1:length(dEs)
                        dE = dEs(dendx);
                        for Trkndx = 1:length(Trks)
                            Trk = Trks(Trkndx);
                            if((sse(sndx,rkndx,r0ndx,dendx,Trkndx) ~= 1e10) && (sse(sndx,rkndx,r0ndx,dendx,Trkndx) ~= 0))
                                % already computed
                                continue;
                            end
                            
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
%                                m_o = zeros(N, 5);
%                                m_o(1, :) = m';
                                for t=2:N
                                    dt = (tdays(t) - tdays(t-1)) * 86400;
                                    Ed2 = (Ed(t)+Ed(t-1))/2;
                                    Ew2 = (Ew(t)+Ew(t-1))/2;

                                    m_new = moisture_model_ext2([1,10,100]',Ed2,Ew2,m,rain(t),dt,1e10,S,rk,r0,Trk);
                                    m_t(t, :) = m_new;
                                    m = m_new;

%                                    m_o(t,:) = moisture_model_ext([1,10,100]',Ed2,Ew2,m_o(t-1,:)',rain(t),dt,1e10);

                                end

                                % compare values
                                ndx = isfinite(fm10);
                                sse(sndx,rkndx,r0ndx,dendx,Trkndx) = sse(sndx,rkndx,r0ndx,dendx,Trkndx) + sum((fm10(ndx)-m_t(ndx,2)).^2);

                                % plot
            %                     subplot(311);
            %                     plot(tdays, m_t(:,2),'r-',tdays,m_o(:,2),'b-',tdays,fm10,'go','markersize',4);
            %                     legend('new', 'old');
            %                     ylabel('fuel moisture [g/100g]');
            %                     subplot(312);
            %                     plot(tdays,rain);
            %                     ylabel('rain intensity [mm/h]');
            %                     subplot(313);
            %                     plot(tdays,Ed,'r-',tdays,Ew,'b-');
            %                     xlabel('time [days in year]');
            %                     pause;
                            end
                            fprintf(' [%g]\n', sse(sndx,rkndx,r0ndx,dendx,Trkndx));

                        end
                    end
                end
            end
        end

        [~,a] = min(sse(:));
        [i1,i2,i3,i4,i5] = ind2sub(size(sse),a);

        fprintf('Minimum reached for : S=%g rk=%g r0=%g dE=%g Trk=%g with error %g\n', Ss(i1),rks(i2),r0s(i3),dEs(i4),Trks(i5),sse(a));

        done = true;
        if(i1==length(Ss)), Ss = [Ss,Ss(end)+0.2]; s = size(sse); s(1) = s(1) + 1; sse2 = zeros(s); sse2(2:end,:,:,:,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i2==1 && rks(1) > 1), rks = [rks(1)-1,rks]; s = size(sse); s(2) = s(2) + 1; sse2 = zeros(s); sse2(:,2:end,:,:,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i2==length(rks)), rks = [rks,rks(end)+1]; s = size(sse); s(2) = s(2) + 1; sse2 = zeros(s); sse2(:,2:end,:,:,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i3==1 && r0s(1) > 0.02), r0s = [r0s(1)-0.02,r0s]; s = size(sse); s(3) = s(3) + 1; sse2 = zeros(s); sse2(:,:,2:end,:,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i3==length(r0s)), r0s = [r0s,r0s(end)+1]; s = size(sse); s(3) = s(3) + 1; sse2 = zeros(s); sse2(:,:,2:end,:,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i4==1), dEs = [dEs(1)-0.01,dEs]; s = size(sse); s(4) = s(4) + 1; sse2 = zeros(s); sse2(:,:,:,2:end,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i4==length(dEs)), dEs = [dEs,dEs(end)+0.01]; s = size(sse); s(4) = s(4) + 1; sse2 = zeros(s); sse2(:,:,:,2:end,:) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i5==1 && Trks(1) > 4), Trks = [Trks(1)-1,Trks]; s = size(sse); s(5) = s(5) + 1; sse2 = zeros(s); sse2(:,:,:,:,2:end) = sse; sse = sse2; done = false; end %#ok<AGROW>
        if(i5==length(Trks)), Trks = [Trks,Trks(end)+1]; s = size(sse); s(5) = s(5) + 1; sse2 = zeros(s); sse2(:,:,:,:,2:end) = sse; sse = sse2; done = false; end %#ok<AGROW>
    end
    
