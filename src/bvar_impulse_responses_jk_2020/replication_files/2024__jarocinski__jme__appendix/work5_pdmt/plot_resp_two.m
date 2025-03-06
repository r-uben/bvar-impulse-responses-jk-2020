function [fh, varminmax] = plot_resp(Cm, Cl, Cu, Cl2, Cu2, vmaturities, vnames, varminmax, b_titles)
% PURPOSE: Plot the impact response of the yield curve and other variables
% INPUTS:
% Cm, Cl, Cu - matrix C (Ns x Nv): median, lower, upper
% vmaturities - Nv vector with maturities of interest rates and NaN for others
vnoty = find(isnan(vmaturities));
vylds = find(~isnan(vmaturities));
xxylds = vmaturities(vylds);
xxnoty = [-0.5 0.5];

xxylds_lab = arrayfun(@(x) sprintf('%dY',x), xxylds, 'UniformOutput', false);
xxylds_lab(logical(rem(xxylds,1))) = {''};
xxylds_lab = [{'0'} xxylds_lab];

%col1 = {[0.7, 0.8, 1], 'EdgeColor', 'none'};
col1 = {[0.5 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none'};
col2 = {[0.5 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 0.4*[1 1 1]}; 


[Ns, Nv] = size(Cm);

if nargin < 8
    scl = 1.01;
    varminmax = [min([Cl*scl; Cu*scl; -ones(1,Nv)]); max([Cl*scl; Cu*scl; ones(1,Nv)])]';
end

if nargin < 9
    b_titles = 1;
end

pos = [5, 1, 6+max(6,3*length(vnoty)), Ns*3.3];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:Ns
    % plot the yield curve
    subplot(Ns, 2, sub2ind([2,Ns], 1, ss))
    toplot = [Cm(ss,vylds); Cl(ss,vylds); Cu(ss,vylds)];
    hold on
    fill([xxylds fliplr(xxylds)], [toplot(2,:) fliplr(toplot(3,:))], col1{:})
    if ~isempty(Cl2)
        %plot(xxylds, Cl2(ss,vylds), ':k', 'LineWidth', 1.25)
        %plot(xxylds, Cu2(ss,vylds), ':k', 'LineWidth', 1.25)
        fill([xxylds fliplr(xxylds)], [Cu2(ss,vylds) fliplr(Cl2(ss,vylds))], col2{:})
    end
    plot(xxylds, toplot(1,:), '-b', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize',14)
    ylim([floor(min(varminmax(vylds,1))) ceil(max(varminmax(vylds,2)))])
    yline(0)
    xticks([0 xxylds])
    xticklabels(xxylds_lab)
    ylabel('basis points')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    if ss==1 && b_titles, title('yield curve'), end
    % plot the other variables
    for ivv = 1:min(2,length(vnoty))
        vv = vnoty(ivv);
        toplot = [Cm(ss,vv); Cl(ss,vv); Cu(ss,vv)];
        subplot(Ns, 4, sub2ind([4,Ns], 2+ivv, ss))
        hold on
        fill([xxnoty fliplr(xxnoty)], [toplot(2) toplot(2) toplot(3) toplot(3)], col1{:})
        if ~isempty(Cl2)
            %plot(xxnoty, Cl2(ss,[vv vv]), ':k', 'LineWidth', 1.25)
            %plot(xxnoty, Cu2(ss,[vv vv]), ':k', 'LineWidth', 1.25)
            fill([xxnoty fliplr(xxnoty)], [Cu2(ss,[vv vv]) Cl2(ss,[vv vv])], col2{:})
        end
        plot(xxnoty, [toplot(1) toplot(1)], '-b', 'LineWidth', 1.5)
        yline(0)
        ylim(varminmax(vv,:))
        xlim([-1 1])
        xticks([])
        grid on
        if ss==1 && b_titles, title(vnames{vv}, 'Interpreter', 'none'), end
    end
end