function suprabasal_plots()

oesophagus_data;

r = 0.135; % ± 0.001
rho = 0.658; % ± 0.005
lambda = 1.32; % ± 0.00 / week
gamma = rho/(1-rho);
mu = rho / 0.85;

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);
margin_left = 0.1; margin_bottom = 0.1;
hpadding = 0.1; vpadding = 0.1;

for i=1:numel(ts3)
    normal{i}(0+1,:) = 0;
    normal{i}(1+1,0+1) = 0;
    
    [m,n] = size(normal{i});
    total = zeros(m+n-1,1);
    for j=1:m
        for k=1:n
            total(j+k-1) = total(j+k-1)+normal{i}(j,k);
        end
    end
    
    raw_theory = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts3(i), z,z,z), ...
        numel(total)-1, 1e-5, 1e-5);
    suprabasal = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts3(i), 0,0,z), ...
        numel(total)-1, 1e-5, 1e-5);
    basal = inverse_z_transform(...
        @(z) generating_function3_shed(r, gamma, mu, lambda*ts3(i), z,z,0), ...
        1, 1e-5, 1e-5);
    Z = 1 - generating_function3_shed(r, gamma, mu, lambda*ts3(i), 0,0,1) - basal(1+1);
    raw_theory = raw_theory - suprabasal;
    raw_theory(1+1) = raw_theory(1+1) - basal(1+1);
    theory = raw_theory / Z;

    i_ = mod(i-1,2); j_ = floor((i-1)/2);
    ah = subplot('Position', [
        margin_left+i_/2*(1-margin_left-hpadding)+i_*hpadding
        margin_bottom+(1-j_)/2*(1-margin_bottom-vpadding)+(1-j_)*vpadding
        (1-margin_left-hpadding)/2 
        (1-margin_bottom-vpadding)/2]);
    set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
    set(ah, 'Box', 'on');
    plot(ah, 2:(m+n-2),total(3:end), 's', 2:(m+n-2),theory(3:end)*sum(total), '-', ...
        0:(m+n-2), sum(total) - binoinv(1-0.158655254, sum(total), 1-theory), '--', ...
        0:(m+n-2), binoinv(1-0.158655254, sum(total), theory), '--');
    xlabel(ah, sprintf('%.1f weeks', ts3(i)), 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ymax = get(ah, 'YLim');
    set(ah, 'YLim', [0 ymax(2)]);
    drawnow;
    
end

print(fh, '-dpdf', '-painters', 'suprabasal-plots');
close(fh);


end