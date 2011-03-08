function r = MDIS(freq, prob)

tot = 0;
for i = 1:numel(freq)
    tot = tot + freq(i);
end

g = freq / tot;

r = tot * kullback_liebler_divergence(g, prob);

end