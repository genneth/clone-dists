function basal

oesophagus_data;

% remove singles
for i=1:numel(ts)
    basal{i}(1+1) = 0;
end

stream = RandStream('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(stream);
sample2(...
    @()(random('beta', 2, 2) / 2), ...
    @()(random('logn', log(1), log(2))), ...
    @()(random('logn', 0, log(2))), ...
    ts2, basal, 100, 100, ...
    'basal.mat');

end