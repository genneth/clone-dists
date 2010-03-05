function oes_compare

t1 = 3/7;
data1 = [1 142; 2 71; 3 8; 4 3];

t2 = 10/7;
data2 = [1 111; 2 109; 3 35; 4 15; 5 4; 6 2];

t3 = 3;
data3 = [1 140; 2 79; 3 38; 4 25; 5 11; 6 4; 7 1; 8 1; 11 1];

t4 = 6;
data4 = [1 73; 2 60; 3 44; 4 26; 5 22; 6 10; 7 3; 8 4; 10 2; 11 4; 13 2; 14 1; 20 1; 28 1];

t5 = 12;
data5 = [1 87; 2 87; 3 50; 4 52; 5 26; 6 24; 7 14; 8 8; 9 9; 10 8; 11 6; 12 7; 13 4; 14 3; 15 3; 16 2; 17 2; 20 2; 22 1; 26 1; 29 1];

t6 = 26;
data6 = [1 40; 2 59; 3 43; 4 18; 5 27; 6 23; 7 16; 8 12; 9 14; 10 10; 11 11; 12 5; 13 8; 14 7; 15 6; 16 3; 17 4; 18 3; 19 3; 20 4; 21 1; 22 2; 23 1; 24 3; 25 4; 27 1; 30 1; 33 1; 34 1; 35 5; 36 1; 37 2; 39 1; 50 2; 56 1; 85 1; 90 1];

t7 = 52;
data7 = [1 12; 2 18; 3 15; 4 14; 5 6; 6 11; 7 8; 8 6; 9 9; 10 9; 11 3; 12 4; 13 4; 14 2; 15 4; 16 4; 17 2; 18 4; 19 2; 20 2; 21 3; 22 4; 23 1; 24 5; 26 2; 27 4; 28 3; 29 1; 30 3; 31 3; 32 3; 33 2; 34 4; 35 2; 36 1; 37 1; 38 1; 40 1; 41 1; 42 1; 43 1; 44 1; 45 1; 47 2; 48 1; 49 1; 54 1; 55 3; 58 1; 60 2; 62 1; 63 1; 70 2; 75 1; 80 1; 85 1; 96 1; 110 1; 130 1; 140 1; 175 1];

ts = [t1 t2 t3 t4 t5 t6 t7];
data = {data1 data2 data3 data4 data5 data6 data7};
av = zeros(size(ts));
for i = 1:numel(ts)
    av(i) = dot(data{i}(:,1), data{i}(:,2) / sum(data{i}(:,2)));
end

lambda = 0.7676;
r = 0.1926;
rho = 0.5164;
gamma = lambda * rho / (1-rho);
tau = rho / (r * lambda);

[p0, ts2] = xi(lambda, r, rho, 100, 0);
loglog(newplot(figure), ts, av, '+', ts2, (1 + (lambda/gamma)*(1-exp(-gamma*ts2))) ./ (1 - p0), '-');
set(gca, 'XLim', [0.1 100]);
xlabel('$t$ / weeks', 'Interpreter', 'latex');
ylabel('$\langle n^\textrm{surv} \rangle$', 'Interpreter', 'latex');
saveas(gcf, 'average-count', 'fig');

%     function scp = smooth(cp)
%         scp = (cp + [1; cp(1:end-1)]) / 2;
%     end

gh = newplot(figure);
colour = 'rgbcmyk';
hold all;
for i = 1:numel(ts);
%     plot(gh, data{i}(:,1) / ts(i), smooth(1 - cumsum(data{i}(:,2)) / sum(data{i}(:,2))), strcat('+:', colour(i)));
    stairs(gh, data{i}(:,1) / ts(i), 1 - cumsum(data{i}(:,2)) / sum(data{i}(:,2)), strcat('-', colour(i)));
    ps = condPb2(exact_pops(lambda, r, lambda * rho / (1-rho), ts(i), floor(10*ts(i))+1));
    plot(gh, (0:floor(10*ts(i))) / ts(i), 1 - cumsum(ps) / sum(ps), strcat('-', colour(i)));
end
plot(gh, data{end}(:,1) / ts(end), exp(-data{end}(:,1) / ts(end) * tau), '-r', 'LineWidth', 2);

hold off;

set(gh, 'XLim', [0 10], 'YLim', [1e-3 1], 'YScale', 'log');
xlabel(gh, '$n/t$ / $\textrm{weeks}^{-1}$', 'Interpreter', 'latex');

legend('', '3 days', ...
    '', '10 days', ...
    '', '3 weeks', ...
    '', '6 weeks', ...
    '', '12 weeks', ...
    '', '26 weeks', ...
    '', '52 weeks', ...
    'limit', ...
    'Location', 'NorthEast');

saveas(gh, 'oes-compare', 'fig');

end