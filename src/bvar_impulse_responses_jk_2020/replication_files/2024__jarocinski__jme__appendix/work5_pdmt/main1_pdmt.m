% Estimate Y = UC, U ~ PDMT(v0,vbar)
% Use data augmentation
% q: Simple Gamma proposal
% priors for W,v0,vbar are proper
clear variables, close all
rng('default');

% if draw_vn: estimate both v0 and vbar
% if ~draw_vn: vbar fixed, store quantities for BF evaluation with bridge sampling
draw_vn = false;
% if draw_vn, the vbar below is the starting point
% if ~draw_vn, the vbar below remains fixed:
vbar = 0.4;

% read data
tab = readtable("data.csv");
fprintf('Data from %s to %s, T=%d\n', tab{1,'date'},tab{end,'date'},size(tab,1))
Y = tab{:,2:end};
ynames = tab.Properties.VariableNames(2:end);
ymaturities = [1/12 2 10 NaN];
[T, N] = size(Y);

% read Student-t results for comparison
try load('baseline_tocompare/maxl.mat'); catch maxl.W=zeros(N,N); maxl.C=zeros(N,N); maxl.C1s=zeros(N,N); end

% options for draw_W_given_QY
draw_W_options = optimoptions(@fminunc, 'Display','off', 'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true, 'HessianFcn','objective');

% prior for W
Wmean = eye(N);
Wprec = 1/(200^2)*eye(N^2);
% starting point
W = chol(cov(Y))\eye(N);

% prior for v0
gmean = 0.1; gstd = 1;
bet0 = gstd^2/gmean; alp0 = gmean/bet0;
draw_v0 = true;
% starting point
v0 = 1;

% prior for vbar->vn
if draw_vn
    gmean = 0.1; gstd = 1;
    betn = gstd^2/gmean; alpn = gmean/betn;
end
% starting point
vn = vbar*ones(1,N);

% starting values for q0, qn
q0 = ones(T,1)*1e-2*v0;
qn = ones(T,N)*1e-2*vbar;

kappa = 3; % scalar in the variance of f

gssettings.ndraws = 2000;
gssettings.burnin = 5000; 50000;
gssettings.saveevery = 5; 500;
gssettings.reportevery = 1e3;

chain.W = nan(N, N, gssettings.ndraws);
chain.Q = nan(T, N, gssettings.ndraws);
chain.v0 = nan(1, gssettings.ndraws);
chain.vbar = nan(1, gssettings.ndraws);
chain.accepted = 0;

chain.vbar_alt = 0.1:0.1:1.2;
chain.r_alt_this = nan(gssettings.ndraws,length(chain.vbar_alt));

nalldraws = gssettings.burnin + gssettings.ndraws*gssettings.saveevery;
fprintf('All draws: %d, burnin: %d, save every: %d, saved draws: %d\n', nalldraws, gssettings.burnin, gssettings.saveevery, gssettings.ndraws)
timing_start = datetime('now');

for draw = 1:nalldraws
    U = Y*W;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw Q

    % draw q0

    % parameters of the proposal gamma
    aa = repmat(0.5*(v0+N), T, 1);
    bb = (U.^2).*repmat((v0+vn).^(-1), T, 1);
    bb = sum(bb,2)+1;
    bb = 2*bb.^(-1);

    % logp
    pl = 0.5*v0 - 1;
    psl = 0.5;
    p1 = -0.5*(1 + sum((U.^2).*repmat((v0+vn).^(-1), T, 1), 2));
    logp = @(x) pl*log(x) + psl*sum(log(x + qn),2) + p1.*x;

    q0 = drawmisgamma_gamma_array(logp, aa, bb, q0, 0);

    % draw qn

    % parameters of the proposal gamma
    A = repmat(0.5*vn + 0.5, T, 1);
    B = (U.^2).*repmat((v0+vn).^(-1), T, 1);
    B = B+1;
    B = 2*B.^(-1);

    % logp
    pl = repmat(0.5*vn - 1, T, 1);
    psl = 0.5;
    p1 = -0.5*(1 + (U.^2).*repmat((v0+vn).^(-1), T, 1));
    logp = @(x) pl.*log(x) + psl*log(q0 + x) + p1.*x;

    qn = drawmisgamma_gamma_array(logp, A, B, qn, 0);

    % draw v0
    if draw_v0
    % logp
    p1 = 0.5*sum(log(q0)) - bet0^(-1) - 0.5*T*log(2);
    pl = alp0 - 1;
    pg = -T;
    psl = -0.5*T;
    psm1 = -0.5*sum((U.^2).*(q0+qn));
    logp = @(x) p1*x + pl*log(x) + pg*gammaln(0.5*x) + psl*sum(log(x + vn')) + sum(psm1'.*(x + vn').^-1);

    % dlogp
    d0 = p1;
    dm1 = pl;
    dpsi0 = 0.5*pg;
    dsm1 = psl;
    dsm2 = -psm1;
    dlogp = @(x) d0 + dm1*x.^(-1) + dpsi0*psi(0.5*x) + dsm1*sum((x+vn').^(-1)) + sum(dsm2'.*(x+vn').^(-2));

    mm = fzero(dlogp, [1e-4 40]);
    ll = -dm1*mm^-2 + 0.5*dpsi0*psi(1,0.5*mm) - dsm1*sum((mm+vn').^-2) + sum(-2*dsm2'.*(mm+vn').^-3);

    [v0, v0logacceptprob] = drawmisgamma_modecurv(logp, mm, ll, v0, 0);
    end
    
    % draw vbar
    if draw_vn
    % logp
    p1 = 0.5*sum(sum(log(qn))) - betn^(-1) - 0.5*T*N*log(2);
    pl = alpn - 1;
    pg = -T*N;
    psl = -0.5*T*N;
    psm1 = -0.5*sum(sum((U.^2).*(q0+qn)));
    logp = @(x) p1*x + pl*log(x) + pg*gammaln(0.5*x) + psl*log(v0 + x) + psm1*(v0 + x).^(-1);

    [vbar, vbarlogacceptprob] = drawmisgamma_modecurvnum(logp, vbar, 0);
    vn = vbar*ones(1,N);
    end

    Q = (qn + q0)./(vn + v0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw W
    [W, accepted] = draw_W_given_QYpriorw(Q, Y, W, Wmean, Wprec, kappa, draw_W_options);
    chain.accepted = chain.accepted + accepted;

    % after burnin 1) fix the order and signs of shocks and 2) reduce kappa
    if draw==gssettings.burnin
        P = normalize1(inv(W), 1:3);
        W = W*P';
        kappa = 1;
    end

    % store the draw
    if draw>gssettings.burnin && ~rem(draw-gssettings.burnin, gssettings.saveevery)
        ii = (draw - gssettings.burnin)/gssettings.saveevery;
        chain.W(:,:,ii) = W;
        chain.Q(:,:,ii) = Q;
        chain.v0(ii) = v0;
        chain.vbar(ii) = vbar;
        
        if ~draw_vn
            % compute and store the quantities for Bayes Factors
            logprior = sum(log(chi2pdf(qn, vbar)),'all');
            logf = loglpdmt(Y, W, q0, qn, v0, vn);

            for jj = 1:length(chain.vbar_alt)
                vbar_jj = chain.vbar_alt(jj);

                logprior_jj = sum(log(chi2pdf(qn, vbar_jj)),'all');
                logf_jj = loglpdmt(Y, W, q0, qn, v0, vbar_jj*ones(1,N));

                chain.r_alt_this(ii,jj) = exp(logprior_jj + logf_jj - logprior - logf);
            end
        end
    end

    if draw==1000 || ~rem(draw, gssettings.reportevery)
        disp(timing_message(draw, nalldraws, timing_start))
        fprintf('Acceptance rate: %.4f\n', chain.accepted/draw)
    end

end
disp(timing_message(draw, nalldraws, timing_start))
if draw_vn
    save(mfilename)
else
    save(mfilename + sprintf("_%0.1f.mat", vbar))
end
fprintf('chain.W(1,1,end)=%f\n', chain.W(1,1,end))


%% report the chain
fprintf('Acceptance rate: %.4f\n', chain.accepted/nalldraws)
fprintf('chain.W(1,1,end)=%f\n', chain.W(1,1,end))

out_path = 'out/';
mkdir(out_path)
pos = [5, 1, 11*N^.6, 8*N^.6];
fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.W(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        plot(x)
        yline(maxl.W(vv, ss))
        if ss==1, ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold'), end
        if vv==1, title(sprintf('u%d',ss), 'FontWeight', 'bold'), end
    end
end
sgtitle('Trace plots of W')
exportgraphics(fh, [out_path 'W_traceplots.pdf'])

fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.W(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        autocorr(x, 'NumLags', min(100,fix(gssettings.ndraws/5)))
        ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold')
        title(sprintf('u%d',ss), 'FontWeight', 'bold')
    end
end
sgtitle('Autocorrelation of draws of W')

fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.W(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        histogram(x)
        xline(maxl.W(vv, ss), 'LineWidth', 1.5)
        if ss==1, ylabel(ynames{vv}, 'Interpreter', 'none', 'FontWeight', 'bold'), end
        if vv==1, title(sprintf('u%d',ss), 'FontWeight', 'bold'), end
    end
end
sgtitle('Histograms of W')
exportgraphics(fh, [out_path 'W_hist.pdf'])

fh = figure;
plot(chain.v0');
sgtitle('Trace plot of v0')
exportgraphics(fh, [out_path 'v0_traceplots.pdf'])

fh = figure;
plot(chain.vbar);
sgtitle('Trace plot of vbar')
exportgraphics(fh, [out_path 'vbar_traceplots.pdf'])

fh = figure;
hold on
h = histogram(chain.v0, 'Normalization', 'pdf');
x = 0:0.01:2;
plot(x, gampdf(x, alp0, bet0))
title('Histogram of v0')
exportgraphics(fh, [out_path 'v0_hist.pdf'])

fh = figure;
histogram(chain.vbar)
sgtitle('Histogram of vbar')
exportgraphics(fh, [out_path 'vbar_hist.pdf'])

fh = figure;
t = 36; n = 2;
plot(squeeze(chain.Q(t,n,:)))
title(sprintf('Trace plot of Q(%d,%d)',t,n))
yline(0)

fh = figure;
plot(log(chain.r_alt_this))
title('Trace plot of log rA/B')

