function normal % normal, suprabasals

oesophagus_data;

% take out unreliable entries
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
    @()(random('beta', 2, 2) / 2), ... % r
    @()(random('logn', log(1), log(2))), ... % gamma
    @()(random('logn', 0, log(2))), ... % lambda
    0.82, ... % suprabasal:basal ratio (m)
    ts3, normal, 1000, ...
    'normal.mat');

end
