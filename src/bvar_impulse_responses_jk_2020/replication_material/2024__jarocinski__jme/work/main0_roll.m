%% Rolling sample maximum likelihood estimation (Figure F.2)
main0_maxlik

Twin = 120;
Nwin = T - Twin + 1;
Crolling = nan(N, N, Nwin);
vrolling = nan(Nwin,N);
exitflags = nan(Nwin,1);
%options.OptimalityTolerance = 1e-9;
%options.FunctionTolerance = 1e-9;
%options.StepTolerance = 1e-09;
options.Display = 'final';
par0 = parmaxlik;

for ww = 1:Nwin
    Ywin = Y(ww:ww+Twin-1,:);
    fun = @(par) nlogl_iidstudent_Wzn(Ywin, reshape(par(1:N^2),N,N), par(N^2+1:end));
    [par,fval,exitflag] = fminunc(fun, par0, options);
    Cwin = inv(reshape(par(1:N^2), N, N));
    Crolling(:,:,ww) = Cwin;
    vrolling(ww,:) = exp(par(N^2+1:end))';
    exitflags(ww) = exitflag;
    par0 = par;
end

tt = tab.date(Twin:end);

for ss = 1:N
    fh = figure;
    x = permute(Crolling(ss,:,:), [3,2,1]);
    if ss==1, x = -x; end
    plot(tt, x, 'LineWidth', 1.5)
    %xline(tt(end-Twin+1))
    %for i = 1:length(tterr), xline(tt2(i)), end
    grid on
    legend(ynames, 'Location', 'Best')
    title(sprintf('u%d',ss))
    exportgraphics(fh, sprintf('%srolling%d-u%d.pdf', out_path, Twin, ss))
end

fh = figure;
plot(tt, vrolling, 'LineWidth', 1.5)
grid on
title('v over rolling samples')
legend("v" + [1:4], "Location", "best")
exportgraphics(fh, sprintf('%srolling%d-v.pdf', out_path, Twin))

