function [xnew, logacceptprob] = drawmisgamma_gamma(logp, aa, bb, x, b_plot)
% PURPOSE: Draw xnew from the distribution given by logp (univariate)
% using the Metropolized Independence Sampler based on Gamma approximation.
% INPUTS:
% logp - function evaluating the log target density
% aa,bb - parameters of the gamma approximation
% x - previous draw
% b_plot - TRUE/FALSe: plot the approximation?
% OUTPUTS:
% xnew - draw of x (might be equal to the previous draw)
% logacceptprob - log of the acceptance probability used for accepting
% 
% Marek Jarocinski 2022-Mar

if b_plot % plot the approximation
    acceptprob = @(xprime) exp(logp(xprime) - logp(x) + log(gampdf(x, aa, bb)) - log(gampdf(xprime, aa, bb)));
    xx = 1e-4:0.01:10;
    figure;
    subplot(2,1,1)
    temp = -max(logp(xx));
    plot(xx, exp(temp+logp(xx)), '-k', 'LineWidth', 2)
    yyaxis right
    plot(xx, gampdf(xx, aa, bb), ':r', 'LineWidth', 2)
    legend('target','approx')
    subplot(2,1,2)
    plot(xx, acceptprob(xx), 'LineWidth', 2);
    xline(x)
    yline(1)
    yline(0)
    ylim([0 1.5])
end

xprime = gamrnd(aa, bb);

logacceptprob = logp(xprime) - logp(x) + log(gampdf(x, aa, bb)) - log(gampdf(xprime, aa, bb));

if log(rand) < logacceptprob
    xnew = xprime;
else
    xnew = x;
end
end