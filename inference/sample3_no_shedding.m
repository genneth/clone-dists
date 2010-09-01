function samples = sample3_no_shedding(rfun, gammafun, lambdafun, ts2, data2, ts3, data3, n)

% lambdafun is unitful; ts2 and ts3 are unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m (Double, Double)
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{prob, r, lambda, gamma, [log p2], [log p3]}]
%     p2 and p3 are likelihoods, not posterior distributions
% n -- number of samples to take

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');
assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');
assert(numel(ts3) > 0, 'no suprabasal data. should run sample2');

% we use the fact that per invocation of clone_dist_* we can extract many
% separate distributions for different lambdas
mult = floor(sqrt(n)/numel(ts3));

wbh = waitbar(0.0, 'initialising...');

samples = cell(n*mult, 6);

maxN2 = 0;
for i = 1:numel(ts2)
    maxN2 = max(maxN2, max(size(data2{i})));
end

maxM3 = 0;
maxN3 = 0;
for i = 1:numel(ts3)
    s = size(data3{i});
    maxM3 = max(maxM3, s(1));
    maxN3 = max(maxN3, s(2));
end
maxK = maxM3 + maxN3;

[Tl Tr Tg P0] = clone_dist_bs_expv_setup(maxK);

tstart = tic;

for i = 1:n
    
    [rp, r] = rfun();
    lambdap = zeros(mult, 1);
    lambda  = zeros(mult, 1);
    for j = 1:mult
        [lambdap(j), lambda(j)] = lambdafun();
    end
    [lambda, SI] = sort(lambda);
    lambdap = lambdap(SI);
    [gammap, gamma] = gammafun();

%     fprintf(1, 'parameters: r = %f, gamma = %f\n', r, gamma);
%     fprintf(1, '            lambdas =\n');
%     fprintf(1, '                %3f\n', lambda);

    lts = zeros(numel(ts3) * mult, 1);
    tsi = zeros(numel(ts3) * mult, 1);
    for j = 1:numel(ts3)
        lts((mult * (j-1) + 1):(mult * j)) = ts3(j) .* lambda;
        tsi((mult * (j-1) + 1):(mult * j)) = j*ones(mult, 1);
    end
    [lts, SI] = sort(lts);
    tsi     = tsi(SI);
    
    dist3_ = clone_dist_bs_expv_cond(Tl, Tr, Tg, P0, r, gamma, lts);
    [m_,n_,~] = size(dist3_);
    dist3 = zeros(m_,n_,mult,numel(ts3));

    for j = 1:numel(ts3)
        ind = find(tsi == j);
        dist3(:,:,:,j) = dist3_(:,:,ind);
    end

    for j = 1:mult
        samples{(i-1)*mult+j,1} = rp*lambdap(j)*gammap;
        samples{(i-1)*mult+j,2} = r; 
        samples{(i-1)*mult+j,3} = gamma;
        samples{(i-1)*mult+j,4} = lambda(j); 

        if numel(ts2) > 0
            samples{(i-1)*mult+j,5} = zeros(numel(ts2), 1);
        else
            samples{(i-1)*mult+j,5} = 1;
        end
        samples{(i-1)*mult+j,6} = zeros(numel(ts3), 1);

        for k = 1:numel(ts3)
            [row,col,v] = find(data3{k});
            for l = numel(v)
                samples{(i-1)*mult+j,6}(k) = ...
                    samples{(i-1)*mult+j,6}(k) ...
                  + log(dist3(row(l),col(l),j,k)) * data3{k}(row(l),col(l));
            end
        end
    end
    
    waitbar(i/n, wbh, sprintf('%.0f completed, %.1f min to go', i/n*100, (toc(tstart)/i*n)/60));

end

close(wbh);

end
