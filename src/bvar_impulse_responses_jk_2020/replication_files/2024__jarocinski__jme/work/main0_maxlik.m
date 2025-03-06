% Estimate Y = UC, Student-t, unknown v
% where
% Y - observable data, i.i.d.
% U - i.i.d. Student-t WITH UNKNOWN v
% C - matrix with the impacts
% We will parameterize in terms of W and z, such that
% YW = U, then recover C = inv(W)
% z = log(v), then recover v = exp(z)
%
% Maximum Likelihood
% Depends: nlogl_iidstudent_Wzn.m (->kr.m), normalize1.m
clear all, close all

% read data
tab = readtable('data.csv');
tab.date.Format = "uuuu-MM-dd";

% Uncomment lines below one-by-one for the shocks reported in Table 1.
%tab = tab(year(tab.date)<2005,:);
%tab = tab(year(tab.date)>2004,:);
%tab(tab.date == datetime('2001-01-03'), :) = [];
%tab(tab.date == datetime('2009-03-18'), :) = [];

fprintf('Data from %s to %s, T=%d\n', tab{1,'date'},tab{end,'date'},size(tab,1))
Y = tab{:,2:end};
ynames = tab.Properties.VariableNames(2:end);
ymaturities = [1/12 2 10 NaN];
[T, N] = size(Y);
nv = N; % 1 or N

out_path = '../results/';
mkdir(out_path)

figure()
for n = 1:N
subplot(2,2,n)
plot(tab.date, Y(:,n))
title(tab.Properties.VariableNames{n+1})
end

% Maximum likelihood / posterior mode estimation of W(->C), z(->v)
options = optimoptions('fminunc');
options.Display = 'iter';
options.MaxFunctionEvaluations = 1e4*N^2;
options.OptimalityTolerance = 1e-8;
options.Algorithm = 'trust-region';
options.SpecifyObjectiveGradient = true;
% starting point
W = inv(chol(cov(Y)) + randn(N,N));
par0 = [W(:); repmat(log(3),N,1)];

% maximum likelihood
fun = @(par) nlogl_iidstudent_Wzn(Y, reshape(par(1:N^2),N,N), par(N^2+1:end));
[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);
W = reshape(parmaxlik(1:N^2), N, N);
C = inv(W);
v = exp(parmaxlik(N^2+1:end));

% Normalize
% 1. Sign: positive interest rate shock
% 2. Order by the maximum impact on consecutive variables
P = normalize1(C, ~isnan(ymaturities))

%P*C
Wnormalized = W*P';
vnormalized = abs(P)*v;


% Minimize again to get the reordered hessian
par0 = [Wnormalized(:); log(vnormalized)];
[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);

maxl.Y = Y;
maxl.ll = -fval;
maxl.parmaxlik = parmaxlik;
maxl.hessian = hessian;
maxl.v = exp(parmaxlik(N^2+1:end));
maxl.W = reshape(parmaxlik(1:N^2), N, N);
maxl.C = inv(maxl.W);
maxl.U = Y*maxl.W;
maxl.U1s = maxl.U./std(maxl.U);
maxl.C1s = diag(std(maxl.U))*maxl.C;
maxl.vdec = (maxl.C1s.^2)./sum(maxl.C1s.^2);
if N==4, temp = [1 2 3 2]; else, temp = 1:N; end
temp = diag(diag(maxl.C(:,temp)));
maxl.C1bp = temp\maxl.C;
maxl.U1bp = maxl.U*temp;

disp('C'), disp(maxl.C)
disp('C1s'), disp(maxl.C1s)
disp('C1bp'), disp(maxl.C1bp)
disp('vdec'), disp(maxl.vdec)
disp('v'), disp(maxl.v)

% save data, maxl, shocks
save(out_path + "data", "tab");
save(out_path + "maxl", "maxl");
u = round(maxl.U, 5);
writetable([tab(:,'date'), array2table(u)], out_path + "U.csv");
u = round(maxl.U1s, 5);
Utable = [tab(:,'date'), array2table(u)];
writetable(Utable, out_path + "U1s.csv");
u = round(maxl.U1bp, 5);
Utable = [tab(:,'date'), array2table(u)];
writetable(Utable, out_path + "U1bp.csv");

% Variance of C, v
JacobianWC = -kron(maxl.W',maxl.W);
Jacobianzv = diag(maxl.v.^-1);
JacobianWzCv = blkdiag(JacobianWC, Jacobianzv);
asyvarCv = inv(JacobianWzCv'*hessian*JacobianWzCv);
asystdCv = sqrt(diag(asyvarCv));
Cstd = reshape(asystdCv(1:N^2),N,N);
vstd = asystdCv(N^2+1:end);
disp(array2table([maxl.v vstd], ...
    "VariableNames", ["MLE","asy-std"], "RowNames", "v"+(1:N)))

Cm = maxl.C;
Cl = Cm - 2*Cstd;
Cu = Cm + 2*Cstd;
fh = plot_resp(diag(std(maxl.U))*Cm, diag(std(maxl.U))*Cl, diag(std(maxl.U))*Cu, ymaturities, ynames);
exportgraphics(fh, out_path + "C1smaxlik_band.pdf")

