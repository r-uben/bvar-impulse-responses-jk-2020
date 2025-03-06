% Estimate Y = UC, Student-t, unknown v
% where
% Y - observable data, i.i.d.
% U - i.i.d. Student-t WITH UNKNOWN SHOCK-SPECIFIC v
% C - matrix with the impacts
% We will actually estimate W in YW = U, then recover C = inv(W)
%
% Metropolis-Hastings sampling
% Depends: nlogl_iidstudent_Wzn or nlogpost_iidstudent_Wzn (->kr.m),
%          timing_message.m
% To be run after main0_maxlik.m and main0_postmode.m
main0_maxlik % loads data, finds maximum likelihood
main0_postmode % finds the posterior mode

if 1
    logtargetdens = @(Y, W, z) -nlogpost_iidstudent_Wzn(Y, W, z, prior);
    maxx = maxp;
else
    logtargetdens = @(Y, W, z) -nlogl_iidstudent_Wzn(Y, W, z);
    maxx = maxl;
end
ndraws = 2000;
saveevery = 5000;
reportevery = 1e6;

% starting point
x = maxl.parmaxlik;
log_dens_x = maxx.ll;

jump_scale = 2/(N-0.5);
chol_inv_H = chol(inv(maxx.hessian))';
jump_chol = jump_scale*chol_inv_H;

chain.logv = nan(ndraws,nv);
chain.W = nan(N, N, ndraws);
chain.log_dens = nan(ndraws, 1);
chain.accepted = 0;
chain.par_mode = maxx.parmaxlik;
chain.max_log_dens = maxx.ll;
nalldraws = ndraws*saveevery;
timing_start = datetime("now");
for draw = 1:nalldraws
        
    % draw candidate
    x_prime = x + jump_chol*randn(N^2+nv, 1);
    
    % decide to accept or not
    log_dens_candidate = logtargetdens(Y, reshape(x_prime(1:N^2),N,N), x_prime(N^2+1:end));

    % truncate
    if any(exp(x_prime(N^2+1:end))<0.45)
        log_dens_candidate = -Inf;
    end

    if log(rand) < log_dens_candidate - log_dens_x
        x = x_prime;
        log_dens_x = log_dens_candidate;
        chain.accepted = chain.accepted + 1;
        
        % keep track of the mode
        if log_dens_x > chain.max_log_dens
            disp(log_dens_x)
            chain.max_log_dens = log_dens_x;
            chain.par_mode = x;
        end
    end
    
    % store the draw
    if ~rem(draw, saveevery)
        chain.logv(draw/saveevery,:) = x(N^2+1:end);
        chain.W(:,:,draw/saveevery) = reshape(x(1:N^2), N, N);
        chain.log_dens(draw/saveevery) = log_dens_x;
    end
    
    if draw==10000 || ~rem(draw, reportevery)
        disp(timing_message(draw, nalldraws, timing_start))
        fprintf('Acceptance rate: %.4f\n', chain.accepted/draw)
    end

end
chain.v = exp(chain.logv);
% update maxl if the MH chain found a higher maximum
if chain.max_log_dens > maxx.ll
    disp('Higher likelihood found, update maxl!')
end
save(out_path + "chain","chain")
%% report the chain
fprintf('Acceptance rate: %.4f\n', chain.accepted/nalldraws)

fh = figure;
plot(chain.log_dens);
yline(maxx.ll);
title('log likelihood')

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
exportgraphics(fh, out_path + "W_traceplots.pdf")

fh = figure('Units','centimeters','Position',pos);
for ss = 1:N
    for vv = 1:N
        x = squeeze(chain.W(vv, ss, :));
        subplot(N, N, sub2ind([N,N], ss, vv))
        autocorr(x, 'NumLags', min(100,fix(ndraws/5)))
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
exportgraphics(fh, out_path + "W_hist.pdf")

fh = figure;
plot(chain.v);
sgtitle('Trace plot of v')
exportgraphics(fh, out_path + "v_traceplots.pdf")

% Histograms of v

fh = figure;
for n = 1:nv
    subplot(ceil(sqrt(nv)),ceil(sqrt(nv)),n)
    histogram(chain.v(:,n))
    title("$v_"+n+"$","Interpreter","latex")
end
exportgraphics(fh, out_path + "v_hist.pdf")

fh = figure;
for n = 1:nv
    subplot(ceil(sqrt(nv)),ceil(sqrt(nv)),n)
    histogram(chain.v(:,n))
    xline(maxl.v(n), 'LineWidth', 1.5)
    title("$v_"+n+"$","Interpreter","latex")
end
exportgraphics(fh, out_path + "v_hist_ml.pdf")

fh = figure;
for n = 1:nv
    subplot(ceil(sqrt(nv)),ceil(sqrt(nv)),n)
    hold on
    histogram(chain.v(:,n), "Normalization","pdf")
    xlimits = xlim;
    xx = linspace(xlimits(1), xlimits(2), 100);
    plot(xx, gampdf(xx, prior.va(n), prior.vb(n)), '-r', 'LineWidth', 2)
    %xline(maxl.v(n), "LineWidth", 1.5)
    title("$v_"+n+"$","Interpreter","latex")
end
exportgraphics(fh, out_path + "v_hist_draws.pdf")

fh = figure;
for n = 1:nv
    subplot(ceil(sqrt(nv)),ceil(sqrt(nv)),n)
    hold on
    histogram(chain.v(:,n), "Normalization","pdf")
    xlimits = xlim;
    xx = linspace(0.1, 15, 100);
    plot(xx, gampdf(xx, prior.va(n), prior.vb(n)), '-r', 'LineWidth', 2)
    title("$v_"+n+"$","Interpreter","latex")
end
exportgraphics(fh, out_path + "v_hist_prior.pdf")
