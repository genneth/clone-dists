function chances = binocompare(experiment, theory)

tot = sum(sum(experiment));
chances = arrayfun(@(s,p) binochance(s,tot,p), full(experiment), full(theory));

end
