function compare(gf, databs, data, rho, r, lambda, t)

k = max(data(:,1) + data(:,2));
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

pbs = condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t, P0)));

hold all;
colour = ' rgbcymk';
for i=2:8;
    errors = sqrt(databs(i,:)) / sum(data(:,3));
    errorbar(gf, 0:numel(databs(i,:))-1, databs(i,:) / sum(data(:,3)), ...
        errors, errors, ...
        strcat('+',colour(i)));
    plot(gf, 0:numel(databs(i,:))-1, pbs(i,1:numel(databs(i,:))), strcat('-',colour(i)));
end
hold off;

xlabel(gf, 'suprabasal', 'Interpreter', 'latex');
xlim = get(gf, 'XLim');
xlim(1) = 0;
set(gf, 'XLim', xlim);
legend('1 basal','', '2','', '3','', '4','', '5','', '6','', '7','');

end
