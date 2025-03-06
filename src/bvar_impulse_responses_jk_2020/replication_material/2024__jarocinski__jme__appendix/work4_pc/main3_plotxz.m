% Plot the responses of variables [X Z] in case Y contained principal components.
% Run this after main3_plots.m

ndraws = size(chain.C1s,3);
chain.CC1s = nan(N, size(Y2XZ,2), ndraws);

for i = 1:ndraws
    % CC1s (the effects of standardized shocks on the variables)
    chain.CC1s(:, :, i) = chain.C1s(:,:,i)*Y2XZ;
end

% Plot: CC1s with uncertainty bands
CC1s_m = maxl.C1s*Y2XZ;
CC1s_l = quantile(chain.CC1s, 0.025, 3);
CC1s_u = quantile(chain.CC1s, 0.975, 3);

fh = plot_resp(CC1s_m, CC1s_l, CC1s_u, xzmaturities, xznames);
exportgraphics(fh, sprintf('%sCC1s_band.pdf', out_path))

