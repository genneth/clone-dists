function intermouse_compare

epsilon = 1e-6;
ts = [21/7 21/7+epsilon 21/7+2*epsilon 30/7 30/7+epsilon 30/7+2*epsilon];
databs{1} = [
    0	13	4	0	0
    12	22	8	2	0
    13	14	7	4	1
    4	2	2	1	0
    2	1	0	0	0
    0	0	0	0	0
    1	1	1	0	0
    0	0	0	0	0
    0	0	0	0	0
    0	0	0	0	1
    0	0	0	1	0
    ];
databs{2} = [
    0	12	13	8	3	3	0	0
    7	11	9	0	2	0	0	0
    10	10	4	6	0	0	0	0
    1	7	5	1	2	1	0	0
    1	3	1	1	1	0	1	0
    0	1	1	1	1	0	1	0
    0	1	2	1	0	1	0	1
    0	1	0	0	0	0	1	0
    0	0	0	0	1	0	0	0
    0	0	1	0	0	0	0	0
    0	0	0	0	0	0	0	0
    0	0	0	0	1	0	0	0
    ];
databs{3} = [
    0	44	22	6	4	1	1
    13	26	7	4	1	0	0
    9	9	6	5	0	1	0
    0	2	3	0	1	0	0
    1	2	2	0	0	0	0
    0	0	0	1	0	0	0
    0	2	0	1	0	0	0
    ];

databs{4} = [
    0	11	6	3	0	1	0	0
    4	15	11	2	5	0	0	0
    4	13	10	5	5	2	2	0
    1	1	1	6	5	0	0	0
    1	1	1	3	5	0	0	0
    0	0	0	3	3	2	0	0
    0	1	1	0	0	1	1	0
    0	0	0	1	1	1	0	2
    0	0	1	1	0	0	0	0
    0	0	0	0	1	0	0	0
    ];
databs{5} = [
    0	11	17	5	2	0	1	0	0
    6	13	16	5	3	0	0	0	0
    7	6	7	4	3	0	0	1	0
    2	6	4	1	2	1	0	0	1
    0	2	5	3	2	1	0	0	0
    1	0	2	5	3	0	1	0	0
    0	0	0	0	0	0	0	0	0
    0	1	0	0	0	1	1	1	1
    0	0	0	0	1	1	1	1	0
    0	0	0	0	0	1	0	1	0
    0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0
    ];
databs{5}(16+1,7+1) = 1;
databs{6} = [
    0	7	3	1	1	1	0
    2	7	5	1	1	1	1
    1	7	6	3	0	0	0
    0	1	1	3	0	1	0
    2	2	3	0	0	2	0
    0	0	0	1	0	0	0
    0	1	0	1	0	1	0
    0	0	0	0	0	0	0
    0	0	0	0	0	1	0
    0	0	0	0	1	0	0
    ];

% take out unreliable entries
for i=1:numel(ts)
    databs{i}(0+1,:)   = 0;
    databs{i}(1+1,0+1) = 0;
end

rs = [0.132 0.216 0.131 0.050 0.062 0.163];
rhos = [0.643 0.553 0.681 0.556 0.572 0.626];
lambdas = [0.83 1.44 0.86 3.11 2.97 1.23];

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);
margin_left = 0.1; margin_bottom = 0.1;
hpadding = 0.1; vpadding = 0.1;
drawnow;

for i=1:numel(ts)
    r = rs(i); rho = rhos(i); lambda = lambdas(i);
    gamma = rho/(1-rho);
    mu = rho / 0.85;
    
    [m,n] = size(databs{i});
    total = zeros(m+n-1,1);
    for j=1:m
        for k=1:n
            total(j+k-1) = total(j+k-1)+databs{i}(j,k);
        end
    end
    
    raw_theory = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts(i), z,z,z), ...
        numel(total)-1, 1e-5, 1e-5);
    suprabasal = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts(i), 0,0,z), ...
        numel(total)-1, 1e-5, 1e-5);
    basal = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts(i), z,z,0), ...
        1, 1e-5, 1e-5);
    Z = 1 - generating_function3_shed(r, gamma, mu, lambda*ts(i), 0,0,1) - basal(1+1);
    raw_theory = raw_theory - suprabasal;
    raw_theory(1+1) = raw_theory(1+1) - basal(1+1);
    theory = raw_theory / Z;

    i_ = mod(i-1,3); j_ = floor((i-1)/3);
    ah = subplot('Position', [
        margin_left+i_/3*(1-margin_left-2*hpadding)+i_*hpadding
        margin_bottom+(1-j_)/2*(1-margin_bottom-vpadding)+(1-j_)*vpadding
        (1-margin_left-2*hpadding)/3 
        (1-margin_bottom-vpadding)/2]);
    set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
    set(ah, 'Box', 'on');
    drawnow;
    plot(ah, 2:(m+n-2),total(3:end), 's', 2:(m+n-2),theory(3:end)*sum(total), '-', ...
        0:(m+n-2), sum(total) - binoinv(1-0.158655254, sum(total), 1-theory), '--', ...
        0:(m+n-2), binoinv(1-0.158655254, sum(total), theory), '--');
    xlabel(ah, sprintf('%.1f weeks', ts(i)), 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ymax = get(ah, 'YLim');
    set(ah, 'YLim', [0 ymax(2)]);
    drawnow;
    
end

print(fh, '-depsc2', '-painters', 'intermouse-suprabasal-plots');
close(fh);

end