function samples = sample3_shed(rfun, gammafun, lambdafun, m, ts2, data2, ts3, data3, nsamples, output_file)

% lambdafun is unitful; ts2 and ts3 are unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m Double
% m :: Double -- average number of suprabasal cells per basal cell
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{r, lambda, gamma, [log p2], [log p3]}]
%     p2 and p3 are likelihoods, not posterior distributions
% nsamples -- number of samples take

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');
assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');
assert(numel(ts3) > 0, 'no suprabasal data. should run sample2');

samples = cell(nsamples, 5);

maxN2 = 0;
for i = 1:numel(ts2)
    maxN2 = max(maxN2, numel(data2{i}) - 1);
end

maxM3 = 0;
maxN3 = 0;
for i = 1:numel(ts3)
    [rows,cols,~] = find(data3{i});
    maxM3 = max(maxM3, max(rows) - 1);
    maxN3 = max(maxN3, max(cols) - 1);
end

for i = 1:nsamples
    
    r = feval(rfun);
    lambda = feval(lambdafun);
    gamma = feval(gammafun);
    rho = gamma / (1 + gamma);
    mu = rho / m; % rho = mu m

    samples{i,1} = r; 
    samples{i,2} = gamma;
    samples{i,3} = lambda; 

    dist3 = clone_dist_bs_shed_cond(r, gamma, mu, lambda * ts3, maxM3, maxN3);
    samples{i,5} = zeros(1,numel(ts3));
    for k = 1:numel(ts3)
        [rows,cols,v] = find(data3{k});
        for l = 1:numel(v)
            samples{i,5}(k) = samples{i,5}(k) + ...
                v(l) * log(dist3(rows(l),cols(l),k));
        end
    end
    
    % and again, with the basal data, if we have it
    if numel(ts2) > 0
        dist2 = clone_dist_b(r, gamma, lambda * ts2, maxN2);
        % condition on b>1
        [m_,~] = size(dist2);
        dist2 = dist2 ./ repmat(1 - dist2(0+1,:) - dist2(1+1,:), m_, 1);
        dist2(0+1,:) = 0; dist2(1+1,:) = 0;

        samples{i,4} = zeros(1,numel(ts2));
        for k = 1:numel(ts2)
            [~,cols,v] = find(data2{k});
            for l = 1:numel(v)
                samples{i,4}(k) = samples{i,4}(k) + ...
                    v(l) * log(dist2(cols(l),k));
            end
        end
    else
        samples{i,4} = [];
    end
    
    partial = samples(1:i,:);
    save(output_file, 'partial');
    
end

end
