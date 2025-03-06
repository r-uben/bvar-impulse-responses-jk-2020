% plot daily LP
% tt, toplot

fh = figure('Units','centimeters','Position',[10 10 7 4]);
hold on
fill([tt' flipud(tt)'], [toplot(:,4)' flipud(toplot(:,5))'], [158,188,218]/255, 'EdgeColor', 'none', 'FaceAlpha',0.5)
fill([tt' flipud(tt)'], [toplot(:,2)' flipud(toplot(:,3))'], [158,188,218]/255, 'EdgeColor', 'none', 'FaceAlpha',0.8)
plot(tt,zeros(size(tt)),'-k');
plot(tt,toplot(:,1)','-k','LineWidth',1.5);
xlabel('horizon h (business days)')
axis tight
grid on


% fname = outdir + varname + "-" + shockname + ".pdf";
% exportgraphics(fh, fname)
