% Normalize the Markov chain of W
% Run this after main1.m

ndraws = size(chain.W,3);
chain.Wnorm = nan(N, N, ndraws);
chain.maxprobs = nan(ndraws,1);

Wref = maxl.W; % reference point for the normalization
Wrefvar = maxl.hessian\eye(size(maxl.hessian,1)); Wrefvar = Wrefvar(1:N^2,1:N^2);
Wrefvar = Wrefvar*1.5;

sss = 1:N; % shocks to permute

allperm = perms(sss);
Nperm = size(allperm,1);
norderswaps = 0;

for i = 1:ndraws
    W = chain.W(:,:,i);
    
    WWW = nan(Nperm,N^2);
    for pp = 1:Nperm
        Wp = W;
        Wp(:,sss) = W(:,allperm(pp,:));
        for n = sss
            temp = Wp\Wref(:,n);
            if temp(n)<0
                Wp(:,n) = -Wp(:,n);
            end
        end
        WWW(pp,:) = Wp(:)';
    end
    probs = mvnpdf(WWW, Wref(:)', Wrefvar);
    probs = probs/sum(probs);
    [maxprob, pp] = max(probs); % pp is the highest probability permutation
    chain.maxprobs(i) = maxprob;
    if 0 % draw permutation proportionally to their probability
        pp = find(cumsum(probs)-rand > 0, 1);
    end
    if pp<Nperm
        norderswaps = norderswaps + 1;
    end
    W = reshape(WWW(pp,:), N, N);

    chain.Wnorm(:, :, i) = W;
end
fprintf('norderswaps = %d\n',norderswaps)
fprintf('min(maxprobs) = %0.3f\n',min(chain.maxprobs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot: trace plots of W normalized
pos = [5, 1, 11*N^.6, 8*N^.6];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.Wnorm(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        plot(x)
        yline(maxl.W(vv, ss))
        if ss==1, ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold'), end
        if vv==1, title(sprintf('u%d',ss), 'FontWeight', 'bold'), end
    end
end
sgtitle('Trace plots of Wnorm')
exportgraphics(fh, out_path + "Wnorm_traceplots.pdf")

% Plot: autocorrelations of W normalized
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.Wnorm(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        autocorr(x, 'NumLags', min(100,fix(ndraws/5)))
        ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold')
        title(sprintf('u%d',ss), 'FontWeight', 'bold')
    end
end
sgtitle('Autocorrelation of draws of Wnorm')

% Plot: histograms of W normalized
pos = [5, 1, 11*N^.6, 8*N^.6];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.Wnorm(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        histogram(x)
        xline(Wref(vv, ss), 'LineWidth', 2)
        if ss==1, ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold'), end
        if vv==1, title(sprintf('u%d',ss), 'FontWeight', 'bold'), end
     end
end
sgtitle('histograms of W normalized')
exportgraphics(fh, out_path + "Wnorm_hist.pdf")
