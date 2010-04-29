function fill_plot(gh, X, Y1, Y2, colour, transparency)

axes(gh);

Y2 = [Y1, fliplr(Y2)];
X2 = [X, fliplr(X)];

fh = fill(X2, Y2, 'b');
set(fh, 'EdgeColor', colour, 'FaceColor', colour, ...
    'EdgeAlpha', transparency, 'FaceAlpha', transparency/2);
ah = get(fh, 'Annotation'); leh = get(ah, 'LegendInformation'); set(leh, 'IconDisplayStyle', 'off');

end
