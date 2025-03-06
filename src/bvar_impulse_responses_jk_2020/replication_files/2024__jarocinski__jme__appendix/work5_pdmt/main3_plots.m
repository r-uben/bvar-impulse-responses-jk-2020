% Process the chain, report the posteriors of quantities derived from W
% Run this after main1 or main2.m
%
% Depends: plot_resp.m

ndraws = size(chain.W,3);
chain.C = nan(N, N, ndraws);
chain.C1s = nan(N, N, ndraws);

for i = 1:ndraws
    if isfield(chain, 'Wnorm')
        W = chain.Wnorm(:,:,i);
    else
        W = chain.W(:,:,i);
    end

    C = inv(W);
    chain.C(:, :, i) = C;   
        
    % C1s (the effects of standardized shocks)
    U = Y*W;
    C1s = diag(std(U))*C;
    chain.C1s(:, :, i) = C1s;
    
end

% Plot: trace plots of C
pos = [5, 1, 11*N^.6, 8*N^.6];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.C(ss, vv, :));
        subplot(N, N, sub2ind([N,N], vv, ss))
        plot(x)
        yline(maxl.C(ss, vv))
        if ss==1, title(ynames{vv}, 'Interpreter', 'none'), end
        if vv==1, ylabel(sprintf('u%d',ss), 'FontWeight', 'bold'), end
     end
end
sgtitle('Trace plots of C')
exportgraphics(fh, out_path + "C_traceplots.pdf")

% Plot: histograms of C
pos = [5, 1, 11*N^.6, 8*N^.6];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.C(ss, vv, :));
        subplot(N, N, sub2ind([N,N], vv, ss))
        histogram(x)
        xline(maxl.C(ss, vv), 'LineWidth', 2)
        if ss==1, title(ynames{vv}, 'Interpreter', 'none'), end
        if vv==1, ylabel(sprintf('u%d',ss), 'FontWeight', 'bold'), end
     end
end
sgtitle('histograms of C')
exportgraphics(fh, out_path + "C_hist.pdf")

% Plot: C with uncertainty bands
C_m = maxl.C;
C_l = quantile(chain.C, 0.025, 3);
C_u = quantile(chain.C, 0.975, 3);
[fh, varminmax] = plot_resp(C_m, C_l, C_u, ymaturities, ynames);
exportgraphics(fh, out_path + "C_band.pdf")

% Plot: C1s with uncertainty bands
C1s_m = maxl.C1s;
C1s_l = quantile(chain.C1s, 0.025, 3);
C1s_u = quantile(chain.C1s, 0.975, 3);
[fh, varminmax] = plot_resp(C1s_m, C1s_l, C1s_u, ymaturities, ynames);
exportgraphics(fh, out_path + "C1s_band.pdf")

% Plot: C1s with uncertainty bands - custom shock reordering
sorder = 1:N;
%sorder = [1 2 3 4]; 
writematrix(sorder, out_path + "sorder.txt");
OO = eye(N); OO = diag(sign(sorder))*OO(abs(sorder),:);
C1s_m1 = OO(1:min(N,4),:)*C1s_m;
C1s_l1 = OO(1:min(N,4),:)*C1s_l; C1s_u1 = OO(1:min(N,4),:)*C1s_u;
[fh, varminmax] = plot_resp(C1s_m1, C1s_l1, C1s_u1, ymaturities, ynames);
exportgraphics(fh, out_path + "C1s_band1.pdf")
if N>4 % plot first 4, then rest
    C1s_m2 = OO(5:end,:)*C1s_m;
    C1s_l2 = OO(5:end,:)*C1s_l; C1s_u2 = OO(5:end,:)*C1s_u;
    fh = plot_resp(C1s_m2, C1s_l2, C1s_u2, ymaturities, ynames, varminmax, 0);
    exportgraphics(fh, out_path + "C1s_band2.pdf")
end

C1s = struct;
C1s.std = OO(1:min(N,4),:)*std(chain.C1s,0,3);
C1s.l = C1s_l1;
C1s.u = C1s_u1;
save(out_path + "C1s.mat", 'C1s')
