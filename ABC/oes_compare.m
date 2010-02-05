function oes_compare

% data from ATRA control
t1 = 3.0+1/7; % 22 days
data1bs = [
	 0  0 39 14  7  4  1  0;
	 0 59 24  6  3  0  0  0;
	32 33 17 15  1  1  0  0;
	 5 11 10  2  3  1  0  0;
	 4  6  3  1  1  0  1  0;
	 5  1  1  2  1  0  1  0;
	 1  4  3  2  0  1  0  1;
	 0  1  0  0  0  0  1  0;
	 0  0  0  0  1  0  0  0;
	 0  0  1  0  1  0  0  0;
	 0  0  0  1  0  0  0  0;
	 0  0  0  0  1  0  0  0
	];
% data1bs(0+1,:) = 0; % condition out suprabasal only
data1 = Pbs2lst(data1bs);

t2 = 3.0+9/7; % 30 days
data2bs = [
	 0  0 26  9  3  2  1  0  0  0;
	 0 35 32  8  9  1  1  0  0  0;
	12 26 23 12  8  2  2  1  0  0;
	 3  8  6 10  7  2  0  0  1  0;
	 3  5  9  6  7  3  0  0  0  1;
	 1  0  2  9  6  2  1  0  0  0;
	 0  2  1  1  0  2  1  0  0  0;
	 0  1  0  1  1  2  1  3  1  0;
	 0  0  1  1  1  2  1  1  0  1;
	 0  0  0  0  2  1  0  1  0  0
	];
% data2bs(0+1,:) = 0; % condition out suprabasal only clones
data2 = Pbs2lst(data2bs);
data2(end+1,:) = [16  7 1];
data2(end+1,:) = [17 10 1];

t = t1;
data = data1;
databs = data1bs;
lambda = 0.78;
r = 0.25;
rho = 0.545;

k = max(data(:,1) + data(:,2));
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

pbs = condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t, P0)));


gf = newplot(figure);
hold all;
colour = 'rgbcymk';
for i=1:7;
    plot(gf, 0:numel(databs(i,:))-1, databs(i,:) / sum(data(:,3)), strcat('+',colour(i)),...
             0:numel(databs(i,:))-1, pbs(i,1:numel(databs(i,:))), strcat('-',colour(i)));
end
hold off;
%set(gf, 'XLim', [0 10], 'YLim', [0 0.3]);

end