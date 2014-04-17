

% This experiment runs all moisture models through one year of data
%


year = 2012;
years = num2str(year);

S = dir(['../data/*fm10*',years,'*']);
N = length(S);

errs = zeros(N,1);
aerrs = zeros(N,1);
sse = zeros(N,1);
counts = zeros(N,1);
snames = cell(N,1);

for i=1:N
    snames{i} = S(i).name(1:5);
    [sd,~,mp] = run_moisture_model_ext(snames{i}, years);
    fm10 = sd.fm10;
    valids = isfinite(fm10);
    errs(i) = sum(fm10(valids) - mp(valids,2));
    aerrs(i) = norm(fm10(valids) - mp(valids,2),1);
    sse(i) = norm(fm10(valids) - mp(valids,2),2);
    counts(i) = sum(valids);
end

save(['station_model_run_',years,'_errors.mat'],'errs','aerrs','sse','counts','snames');
