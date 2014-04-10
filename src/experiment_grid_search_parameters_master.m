

% Master process for experiment_grid_search_parameters_slave
% Constructs a matrix containing all combinations of parameters
% and passes segments to slaves.

% year
year = '2013';

% all combinations
dEs = -0.06:0.01:0.06;
Ss = 0.2:0.1:1.6;
rks = 5:16;
r0s = 0.02:0.01:0.14;
Trks = 6:1:16;

% number of processes
N = 16;


% no of combinations
n = prod([length(dEs),length(Ss),length(rks),length(r0s),length(Trks)]);


items = zeros(n,5);
fprintf('Have %d parameter combinations.\n', size(items,1));

Ni = floor(n / N);

% assign all combinations into items matrix
j = 1;
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
                    items(j,:) = [S,rk,r0,dE,Trk];
                    j = j + 1;
                end
            end
        end
    end
end

save grid_search_items.mat items

% write into a script file matlab invocations
f = fopen('grid_search.sh', 'w');
for i=1:N
    from = (i-1)*Ni+1;
    if(i<N),to=i*Ni; else to=n; end
    fprintf(f, 'nohup matlab -r "from = %d; to = %d; sse = experiment_grid_search_parameters_slave(from,to,%s); save grid_search_results_%d.mat sse from to;" &> grid_search_%d.log &\n', from, to, year, i, i);
end
fclose(f);

