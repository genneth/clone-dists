function single_cell_plot

fh = figure;
ah = newplot(fh);

ts = linspace(0,60,100);

% parameters from clonal labelling
rhoS = 0.044; lambdaS = 1/30;
lambda = 0.7676;
sc = rhoS * exp(-lambdaS * ts) + (1-rhoS) * exp(-lambda * ts);

% 6 layers
sc6 = single_cell(0.2, 1.0, 1/12, lambda * ts);

% 3 layers
sc3 = single_cell(0.2, 1.0, 1/6, lambda * ts);

% 2 layers
sc2 = single_cell(0.2, 1.0, 1/4, lambda * ts);

% 1 layer
sc1 = single_cell(0.2, 1.0, 1/2, lambda * ts);

semilogy(ah, ts, sc, ts, sc6, ts, sc3, ts, sc2, ts, sc1);

legend(ah, {'clonal labelling exp.', '6 layers', '3 layers', '2 layers', '1 layer'});
set(ah, 'FontName', 'Times', 'FontSize', 8);
xlabel(ah, 'time (weeks)', 'FontName', 'Times', 'FontSize', 8);
ylabel(ah, 'proportion', 'FontName', 'Times', 'FontSize', 8);

end
