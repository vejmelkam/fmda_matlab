%
% Run a grid search through the 5d parameter space of [dE,S,rk,r0,Trk].
%
% synopsis: [sse,Nsse] = experiment_grid_search_parameters_wrf(wrf_file,stop_at_da)
%
%

function [sse,Nsse] = experiment_grid_search_parameters_wrf(wrf_file,stop_at_da)

    % read out data from wrf file
    wrf_t = ncread(wrf_file, 'Times')';
    Nt = size(wrf_t,1);
    psfc = ncread(wrf_file,'PSFC');
    q2 = ncread(wrf_file, 'Q2');
    t2 = ncread(wrf_file,'T2');
    glon = ncread(wrf_file,'XLONG');
    glon = glon(:,:,1);
    glat = ncread(wrf_file,'XLAT');
    glat = glat(:,:,1);
    rain = cat(3,zeros(size(glon)),diff(ncread(wrf_file,'RAINC'),1,3) + diff(ncread(wrf_file,'RAINNC'),1,3));
    [ed,ew,relh] = equilibrium_moisture(psfc,q2,t2);
    relh = relh / 100;
    
    assert(all(ed(:) > 0));
    assert(all(ew(:) > 0));
    
    if(any(relh(:) < 0))
        warning('fixing WRF relative humidity, found %d negative values', sum(relh(:) < 0));
        relh(relh < 0) = 0;
    end
    if(any(relh(:) > 1))
        warning('fixing WRF relative humidity, found %d negative values', sum(relh(:) > 1));
        relh(relh > 1) = 1;
    end
    
    % model params kept fixed
    Tk = [1,10,100]';

    % initial grid search
    dEs = -0.1:0.01:0.02;
    Ss = 0.2:0.2:2.0;
    rks = 1:2:14;
    r0s = 0.04:0.02:0.1;
    Trks = 4:2:14;    
    
    % initial grid search
%     dEs = -0.05:0.01:-0.03;
%     Ss = 0.2:0.2:0.6;
%     rks = 1:3;
%     r0s = 0.06:0.02:0.1;
%     Trks = 5:7;
    
    % outer limits beyond which we don't want the optimizer to move
%     dEsl = [-0.1,0.1];
%     Ssl = [0.2,2.6];
%     rksl = [1,16];
%     r0sl = [0.01,0.1];
%     Trksl = [4,18];
    
%     dEsl = [-0.1,0.1];
%     Ssl = [0.2,0.6];
%     rksl = [1,3];
%     r0sl = [0.06,0.08];
%     Trksl = [5,7];
    
    % parse the dates
    Ts = zeros(Nt,1);
    for i=1:Nt
        ti = wrf_t(i,:);
        Ts(i) = datenum(str2double(ti(1:4)),str2double(ti(6:7)),str2double(ti(9:10)),...
                        str2double(ti(12:13)),str2double(ti(15:16)),str2double(ti(18:19)));
        % convert rain to mm/h units
        if(i>1)
            rain(:,:,i) = rain(:,:,i) / ((Ts(i) - Ts(i-1)) * 24);
        end
    end

    % remove any rain that is too light to affect the moisture model
    rain(rain < 0.01) = 0;
    
    % load data from correct year
    year = str2double(ti(1:4));
    years = num2str(year);
    S = dir(['../data/*fm10*',years,'*']);
    N = length(S);
    
    % load all station data & find closest grid point
    Nst = 0;
    have_rain = false;
    sds = cell(N,1);
    for i=1:N
        sds_i = load_station_data_for_wrf(S(i).name(1:5),Ts);
        sds_i.stid = S(i).name(1:5);
        [lati,loni,disti] = find_closest_point(sds_i.lat,sds_i.lon,glat,glon);
        if(disti < 3)
            Nst = Nst + 1;
            sds_i.glat = lati;
            sds_i.glon = loni;
            sds{Nst} = sds_i;
            if(any(rain(loni,lati,:) > 0))
                have_rain = true;
            end
        else
            fprintf('Skipping station %s, too far from grid (%g km)\n', S(i).name(1:5), disti);
        end
    end
    sds = sds(1:Nst);    
    
    if(~have_rain)
        warning('WRF forecast has no rain');
    end
    
    sse = zeros(length(Ss),length(rks),length(r0s),length(dEs),length(Trks));
    Nsse = zeros(size(sse));
    
    done = false;
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
                            
                            % initialize models
                            k=3;                            
                            ms = zeros(Nst,k+2);
                            for s=1:Nst
                                sdi = sds{s};
                                i = sdi.glon;
                                j = sdi.glat;
                                ms(s,1:3) = 0.5 * (ed(i,j,1) + ew(i,j,1));
                                ms(s,4) = dE;
                            end
                                
                            % run models
                            for t=2:stop_at_da
                                for s=1:Nst
                                    sdi = sds{s};
                                    i = sdi.glon;
                                    j = sdi.glat;
                                    Ed2 = 0.5 * (ed(i,j,t) + ed(i,j,t-1));
                                    Ew2 = 0.5 * (ew(i,j,t) + ed(i,j,t-1));
                                    r = rain(i,j,t);
                                    dt = (Ts(t)-Ts(t-1))*86400;
                                    ms(s,:) = moisture_model_ext2(Tk,Ed2,Ew2,ms(s,:)',r,dt,1e10,S,rk,r0,Trk);
                                end
                                
                                % accumulate errors
                                [~,fm10o,~] = find_valid_obs(Ts(t),sds);
                                valids = find(isfinite(fm10o));
                                sse(sndx,rkndx,r0ndx,dendx,Trkndx) = sse(sndx,rkndx,r0ndx,dendx,Trkndx) + sum((ms(valids,2)-fm10o(valids)).^2);
                                Nsse(sndx,rkndx,r0ndx,dendx,Trkndx) = Nsse(sndx,rkndx,r0ndx,dendx,Trkndx) + length(valids);
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

        % if any of the minima is reached at a boundary value, increase the
        % boundary that way and continue
        done = true;
%         s = size(sse);
%         sse2 = zeros(s);
%         Nsse2 = zeros(s);
%         
%         if((i1==length(Ss)) && (Ss(end)<Ssl(2))), Ss = [Ss(2:end),Ss(end)+0.2]; sse2(1:end-1,:,:,:,:) = sse(2:end,:,:,:,:); Nsse2(1:end-1,:,:,:,:) = Nsse(2:end,:,:,:,:); sse = sse2; Nsse = Nsse2; done = false; end
%         
%         if((i2==1) && (rks(1) > rksl(1))), rks = [rks(1)-1,rks(1:end-1)]; sse2(:,2:end,:,:,:) = sse(:,1:end-1,:,:,:); sse = sse2; Nsse2(:,2:end,:,:,:) = Nsse(:,1:end-1,:,:,:); Nsse = Nsse2; done = false; end
%         if((i2==length(rks)) && (rks(end) < rksl(2))), rks = [rks(2:end),rks(end)+1]; sse2(:,1:end-1,:,:,:) = sse(:,2:end,:,:,:); sse = sse2; Nsse2(:,1:end-1,:,:,:) = Nsse(:,2:end,:,:,:); Nsse = Nsse2; done = false; end
%         
%         if((i3==1) && (r0s(1) > r0sl(1))), r0s = [r0s(1)-0.02,r0s(1:end-1)]; sse2(:,:,2:end,:,:) = sse(:,:,1:end-1,:,:); sse = sse2; Nsse2(:,:,2:end,:,:) = Nsse(:,:,1:end-1,:,:); Nsse = Nsse2; done = false; end
%         if((i3==length(r0s)) && (r0s(end) < r0sl(2))), r0s = [r0s(2:end),r0s(end)+0.02]; sse2(:,:,1:end-1,:,:) = sse(:,:,2:end,:,:); sse = sse2; Nsse2(:,:,1:end-1,:,:) = Nsse(:,:,2:end,:,:); Nsse = Nsse2; done = false; end
%         
%         if((i4==1) && (dEs(1) > dEsl(1))), dEs = [dEs(1)-0.01,dEs(1:end-1)]; sse2(:,:,:,2:end,:) = sse(:,:,:,1:end-1,:); sse = sse2; Nsse2(:,:,:,2:end,:) = Nsse(:,:,:,1:end-1,:); Nsse = Nsse2; done = false; end
%         if((i4==length(dEs)) && (dEs(2) < dEsl(2))), dEs = [dEs(2:end),dEs(end)+0.01]; sse2(:,:,:,1:end-1,:) = sse(:,:,:,2:end,:); sse = sse2; Nsse2(:,:,:,1:end-1,:) = Nsse(:,:,:,2:end,:); Nsse = Nsse2; done = false; end
%         
%         if((i5==1) && (Trks(1) > Trksl(1))), Trks = [Trks(1)-1,Trks(1:end-1)]; sse2(:,:,:,:,2:end) = sse(:,:,:,:,1:end-1); sse = sse2; Nsse2(:,:,:,:,2:end) = Nsse(:,:,:,:,1:end-1); Nsse = Nsse2; done = false; end
%         if((i5==length(Trks)) && (Trks(end) < Trksl(2))), Trks = [Trks(2:end),Trks(end)+1]; sse2(:,:,:,:,1:end-1) = sse(:,:,:,:,2:end); sse = sse2; Nsse2(:,:,:,:,1:end-1) = Nsse(:,:,:,:,2:end); Nsse = Nsse2; done = false; end
        
    end
    
end
    
    