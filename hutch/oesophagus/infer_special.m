function infer_special

oesophagus_data;

% take out unreliable entries
for i=1:numel(ts2)
    basal{i}(1+1) = 0;
end

for i=1:numel(ts3)
    normal{i}(0+1,:)   = 0;
    normal{i}(1+1,0+1) = 0;
end

% cut down on the data to make testing faster
% ts2 = ts2(4); basal = basal(4);
% ts3 = ts3(1:2); normal = normal(1:2);

stream = RandStream('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(stream);
sample3_shed(...
    @()(0.25), ...
    @()(random('logn', log(1), log(2))), ...
    @()(random('logn', 0, log(2))), ...
    0.85, ...
    ts2, basal, ts3, normal, 100, ...
    'special.mat');

end
