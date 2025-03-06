% Bayes factor using bridge sampler
% BF = m_A(Y)/m_B(Y)
% model A: 
vbarA = 0.6;
% model B:
vbarB = 0.4;

load(sprintf("main1_pdmt_%0.1f.mat", vbarB));
chainB = chain;
chainB.vbar(1)
jj = chainB.vbar_alt == vbarA;
rAB = chainB.r_alt_this(:,jj);

load(sprintf("main1_pdmt_%0.1f.mat", vbarA));
chainA = chain;
chainA.vbar(1)
jj = chainA.vbar_alt == vbarB;
rBA = chainA.r_alt_this(:,jj);

figure;
subplot(1,2,1)
plot(log(rBA))
title("rB/A")
subplot(1,2,2)
plot(log(rAB))
title("rA/B")

niter = 100;
BFs = nan(niter,1);
BFs(1) = 1;
for i = 2:niter
    numer = mean(rAB./(rAB + BFs(i-1)));
    denom = mean(rBA./(rBA*BFs(i-1) + 1));
    BFs(i) = numer/denom;
end

figure;
plot(log(BFs))

% in this case it is easy to see where the log is converging
logbf = mean(log(BFs(end-1:end)))
bf = exp(logbf)
