function basal_check

oesophagus_data;

% use the suprabasal counts to estimate tau

basal = {};
for i=1:numel(ts3)
    basal{i} = sum(normal{i},2);
    basal{i}(0+1) = 0; % remove extinct ones, because of fake data
    basal{i}(1+1) = 0; % remove singles
end

stream = RandStream('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(stream);
sample2(...
    @()(random('beta', 2, 2) / 2), ...
    @()(random('logn', log(1), log(2))), ...
    @()(random('logn', 0, log(2))), ...
    ts3, basal, 100, 100, ...
    'basal-check.mat');

end